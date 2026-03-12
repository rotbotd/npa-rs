# Lecture 8: Implementing It All

Time to put everything together. This lecture walks through the actual Rust implementation.

## Architecture Overview

```
src/
├── semiring.rs          # Semiring, Admissible traits
├── boolean_matrix.rs    # Dense N×N Boolean matrices
├── sparse_matrix.rs     # Sparse Boolean matrices (for scale)
├── lifted_matrix.rs     # Matrices with explicit ⊥ element
├── expr.rs              # Regular expressions over semirings
├── cfg.rs               # Control flow graphs and dominator trees
├── tarjan.rs            # Path expression algorithm
├── differentiate.rs     # Formal differentiation
├── regularize.rs        # τ_Reg transformation (LCFL → left-linear)
├── npa.rs               # Main NPA-TP solver loop
├── linear.rs            # Fast path for left-linear equations
└── fuzz.rs              # Differential fuzzing against Kleene
```

## The Semiring Trait

```rust
pub trait Semiring: Clone + Eq + std::fmt::Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn combine(&self, other: &Self) -> Self;  // ⊕
    fn extend(&self, other: &Self) -> Self;   // ⊗
    fn star(&self) -> Self;                   // Kleene star
    
    // For sized elements (matrices)
    fn zero_sized(size: usize) -> Self { Self::zero() }
    fn one_sized(size: usize) -> Self { Self::one() }
    fn size(&self) -> usize { 0 }
}
```

The `_sized` variants handle matrices where we need to know dimensions. A 2×2 zero matrix is different from a 4×4 zero matrix.

## The Admissible Trait

```rust
pub trait Admissible: Semiring {
    type Tensor: Semiring;

    fn transpose(&self) -> Self;
    fn tensor(&self, other: &Self) -> Self::Tensor;
    fn detensor_transpose(t: &Self::Tensor) -> Self;
}

pub fn coupling<S: Admissible>(a: &S, b: &S) -> S::Tensor {
    a.transpose().tensor(b)
}
```

The associated type `Tensor` lets us express that tensoring two N×N matrices gives an N²×N² matrix.

## Boolean Matrices

The core domain for predicate abstraction:

```rust
pub struct BoolMatrix {
    n: usize,
    data: Vec<bool>,  // row-major n×n
}

impl Semiring for BoolMatrix {
    fn zero() -> Self { /* all false */ }
    fn one() -> Self { /* identity matrix */ }
    
    fn combine(&self, other: &Self) -> Self {
        // Pointwise OR
        BoolMatrix {
            n: self.n,
            data: self.data.iter()
                .zip(&other.data)
                .map(|(a, b)| *a || *b)
                .collect(),
        }
    }
    
    fn extend(&self, other: &Self) -> Self {
        // Matrix multiply with OR/AND
        let mut result = BoolMatrix::zero_sized(self.n);
        for i in 0..self.n {
            for j in 0..self.n {
                for k in 0..self.n {
                    if self.get(i, k) && other.get(k, j) {
                        result.set(i, j, true);
                    }
                }
            }
        }
        result
    }
    
    fn star(&self) -> Self {
        // Transitive closure via repeated squaring
        // a* = I ⊕ a ⊕ a² ⊕ a⁴ ⊕ ... until stable
        let mut result = BoolMatrix::one_sized(self.n);
        let mut power = self.clone();
        loop {
            let new_result = result.combine(&power);
            if new_result == result {
                return result;
            }
            result = new_result;
            power = power.extend(&power);
        }
    }
}
```

### Transpose

```rust
fn transpose(&self) -> Self {
    let mut result = BoolMatrix::zero_sized(self.n);
    for i in 0..self.n {
        for j in 0..self.n {
            result.set(j, i, self.get(i, j));
        }
    }
    result
}
```

### Kronecker Product (Tensor)

```rust
fn tensor(&self, other: &Self) -> BoolMatrix {
    let n2 = self.n * self.n;
    let mut result = BoolMatrix::zero_sized(n2);
    
    for i1 in 0..self.n {
        for j1 in 0..self.n {
            for i2 in 0..other.n {
                for j2 in 0..other.n {
                    let row = i1 * other.n + i2;
                    let col = j1 * other.n + j2;
                    result.set(row, col, self.get(i1, j1) && other.get(i2, j2));
                }
            }
        }
    }
    result
}
```

### Detensor-Transpose

This is the magic operation. From equation 46 in the paper:

```
(t,·)(T)[A, A'] = ∃B, B' : T[(A',B), (A,B')] ∧ (B = B')
```

In code:

```rust
fn detensor_transpose(t: &BoolMatrix) -> Self {
    let n2 = t.n;
    let n = (n2 as f64).sqrt() as usize;
    
    let mut result = BoolMatrix::zero_sized(n);
    
    for a in 0..n {
        for a_prime in 0..n {
            // Check if exists b such that T[(a', b), (a, b)] is true
            for b in 0..n {
                let row = a_prime * n + b;  // (A', B)
                let col = a * n + b;        // (A, B') where B' = B
                if t.get(row, col) {
                    result.set(a, a_prime, true);
                    break;
                }
            }
        }
    }
    result
}
```

The diagonal constraint `B = B'` is what prevents cross-terms!

## Expressions

Regular expressions over a semiring:

```rust
pub enum Expr<S: Semiring> {
    Zero,
    One,
    Const(S),
    Var(usize),
    Combine(Arc<Expr<S>>, Arc<Expr<S>>),
    Extend(Arc<Expr<S>>, Arc<Expr<S>>),
    Star(Arc<Expr<S>>),
}
```

We use `Arc` for subexpression sharing — identical subtrees share memory.

### Evaluation

```rust
impl<S: Semiring> Expr<S> {
    pub fn eval(&self, vars: &[S]) -> S {
        match self {
            Expr::Zero => S::zero(),
            Expr::One => S::one(),
            Expr::Const(s) => s.clone(),
            Expr::Var(i) => vars[*i].clone(),
            Expr::Combine(a, b) => a.eval(vars).combine(&b.eval(vars)),
            Expr::Extend(a, b) => a.eval(vars).extend(&b.eval(vars)),
            Expr::Star(a) => a.eval(vars).star(),
        }
    }
}
```

## Differentiation

From Lecture 6, the derivative rules:

```rust
pub fn differentiate<S: Semiring>(expr: &Expr<S>, var_idx: usize) -> Expr<S> {
    match expr {
        Expr::Zero | Expr::One | Expr::Const(_) => Expr::Zero,
        
        Expr::Var(i) => {
            if *i == var_idx { Expr::One } else { Expr::Zero }
        }
        
        Expr::Combine(a, b) => {
            differentiate(a, var_idx).combine(differentiate(b, var_idx))
        }
        
        Expr::Extend(a, b) => {
            // Product rule: D[a⊗b] = D[a]⊗b ⊕ a⊗D[b]
            let da = differentiate(a, var_idx);
            let db = differentiate(b, var_idx);
            da.extend((**b).clone()).combine((**a).clone().extend(db))
        }
        
        Expr::Star(a) => {
            // D[a*] = a* ⊗ D[a] ⊗ a*
            let da = differentiate(a, var_idx);
            let a_star = Expr::Star(a.clone());
            a_star.clone().extend(da).extend(a_star)
        }
    }
}
```

The star rule `D[a*] = a* ⊗ D[a] ⊗ a*` looks like it creates sandwiches — but the `a*` terms are **constants** (the current ν), not variables. So the result is still linear in the new variable.

## The τ_Reg Transformation

Converting LCFL to left-linear over the tensor semiring:

```rust
// LCFL term: a ⊗ Y_i ⊗ b
pub struct LcflTerm<S: Semiring> {
    pub left: Expr<S>,   // a
    pub var_idx: usize,  // i (which Y variable)
    pub right: Expr<S>,  // b
}

// Convert to tensor coefficient: a^t ⊗ b
pub fn tau_reg_term<S: Admissible>(term: &LcflTerm<S>, nu: &[S]) -> S::Tensor {
    let a_val = term.left.eval(nu);
    let b_val = term.right.eval(nu);
    coupling(&a_val, &b_val)
}
```

The key transformation: `a ⊗ Y ⊗ b` becomes `Z ⊗_T (a^t ⊗ b)`.

## The NPA-TP Solver

Putting it all together:

```rust
pub struct NpaSolver<S: Admissible> {
    pub n: usize,                              // number of equations
    pub rhs: Vec<Expr<S>>,                     // right-hand sides
    pub dep_graph: Cfg<DepLabel>,              // dependence graph
    pub path_exprs: HashMap<usize, Expr<DepLabel>>,  // from Tarjan
    pub lcfl_terms: Vec<Vec<Vec<LcflTerm<S>>>>,      // precomputed
    pub one: S,                                // identity element
}

impl<S: Admissible> NpaSolver<S> {
    pub fn solve(&self, max_rounds: usize) -> NpaResult<S> {
        // Initialize: ν = f(0)
        let zero_values: Vec<S> = vec![S::zero_sized(self.one.size()); self.n];
        let mut nu: Vec<S> = self.rhs.iter()
            .map(|rhs| rhs.eval(&zero_values))
            .collect();
        
        let mut round = 0;
        
        loop {
            round += 1;
            if round > max_rounds { break; }
            
            // Compute tensor coefficients T_kj at current ν
            let coeffs = self.compute_coefficients(&nu);
            
            // Solve using path expressions from Tarjan
            let z_values = self.solve_with_path_exprs(&nu, &coeffs);
            
            // Apply detensor-transpose to get new ν
            let new_nu: Vec<S> = z_values.iter()
                .map(|z| S::detensor_transpose(z))
                .collect();
            
            // Check convergence
            if new_nu == nu {
                return NpaResult { values: new_nu, rounds: round };
            }
            
            nu = new_nu;
        }
        
        NpaResult { values: nu, rounds: round }
    }
}
```

### Key Insight: Structure is Precomputed

The `path_exprs` from Tarjan don't change between rounds. Only the **leaf values** change. This is huge for efficiency — we're not rebuilding path expressions every round.

### Memoization

When evaluating path expressions, identical subexpressions should return the same value. We use pointer-based memoization:

```rust
fn eval_path_expr_memo<S: Admissible>(
    expr: &Expr<DepLabel>,
    label_values: &HashMap<DepLabel, S::Tensor>,
    memo: &mut HashMap<*const Expr<DepLabel>, S::Tensor>,
) -> S::Tensor {
    let ptr = expr as *const Expr<DepLabel>;
    if let Some(cached) = memo.get(&ptr) {
        return cached.clone();
    }
    // ... compute and cache
}
```

Without memoization, DAG expressions become exponential.

## Sparse Matrices

For real programs, N = 2^|predicates| gets big. We use sparse representation:

```rust
pub struct SparseBoolMatrix {
    pub n: usize,
    pub entries: HashSet<(usize, usize)>,  // only store true entries
}
```

Operations become:
- `combine`: union of entry sets
- `extend`: iterate over pairs that match
- Memory: O(|entries|) instead of O(n²)

## Linear Fast Path

Many dataflow problems are actually left-linear (no sandwiches). For these, we skip the tensor machinery:

```rust
pub fn is_left_linear<S: Semiring>(expr: &Expr<S>) -> bool {
    // Check if variables only appear on the left of extends
    // i.e., always "X ⊗ const" never "const ⊗ X ⊗ const"
}

pub fn solve_linear<S: Semiring>(rhs: Vec<Expr<S>>) -> Vec<S> {
    // Extract coefficient matrix A and constants c
    // Solve X = c ⊕ A·X directly with Kleene star on A
}
```

This covers most intraprocedural dataflow — huge speedup.

## Testing

We fuzz against naive Kleene iteration:

```rust
fn fuzz_against_kleene<S: Admissible + Random>() {
    for _ in 0..1000 {
        let rhs = random_equations();
        let npa_result = solve_npa(rhs.clone());
        let kleene_result = kleene_iteration(rhs);
        assert_eq!(npa_result.values, kleene_result.values);
    }
}
```

If NPA-TP and Kleene agree on random inputs, we're probably correct.

## Performance Considerations

1. **Subexpression sharing**: Use `Arc` and don't duplicate
2. **Memoization**: Cache evaluation results by pointer
3. **Sparse matrices**: Essential for |predicates| > 10
4. **Linear fast path**: Skip tensor for most intraprocedural
5. **Preprocessing**: Tarjan runs once, reuse structure

## Problem Set 8

### Problem 8.1: Implement Boolean Semiring

Write a `BoolSemiring` that implements `Semiring` for plain booleans:
- `zero() = false`
- `one() = true`
- `combine = ||`
- `extend = &&`
- `star() = true`

Test that it satisfies the semiring axioms.

### Problem 8.2: Matrix Multiplication

Implement `extend` for 2×2 Boolean matrices by hand. Given:
```
A = [[1,0],[1,1]]
B = [[0,1],[0,0]]
```
Compute A ⊗ B using OR/AND.

### Problem 8.3: Transitive Closure

For the matrix A = [[0,1,0],[0,0,1],[0,0,0]], compute A* step by step.

How many squaring iterations does it take to converge?

### Problem 8.4: Kronecker Product

Compute the Kronecker product of:
```
A = [[1,0],[0,1]]  (2×2 identity)
B = [[1,1],[0,0]]
```

The result should be 4×4.

### Problem 8.5: Detensor-Transpose

Given a 4×4 tensor matrix T where only T[(0,0), (1,0)] = true:
- What is (t,·)(T)?
- What does this represent in terms of relations?

### Problem 8.6: Expression Simplification

The expression `Zero.combine(x)` should simplify to `x`. The expression `One.extend(x)` should simplify to `x`.

Implement a `simplify` function for `Expr<S>` that applies these identities.

### Problem 8.7: DAG Expressions

Consider the expression built as:
```rust
let x = Expr::var(0);
let y = x.clone().extend(x.clone());  // x²
let z = y.clone().combine(y.clone()); // x² ⊕ x²
```

Without memoization, evaluating `z` would evaluate `x` how many times?

With pointer-based memoization, how many times?

### Problem 8.8: Running the Code

Clone the npa-rs repo and run the tests:
```bash
git clone https://github.com/rotbotd/npa-rs
cd npa-rs
cargo test
```

Pick a test case and trace through the NPA-TP algorithm by adding print statements. Watch Newton converge in ~3 rounds.

---

## What's Next?

You've now seen the full NPA-TP pipeline:

1. **Semirings** abstract the dataflow domain
2. **Path expressions** encode all paths symbolically
3. **Tarjan** computes path expressions efficiently
4. **Newton** linearizes nonlinear fixed-point equations
5. **Tensor products** regularize LCFL systems
6. **Detensor-transpose** recovers the solution

The implementation in npa-rs is a working prototype. For production use, you'd want:
- Better sparse matrix representations (BDDs, tries)
- Incremental updates when programs change
- Integration with actual compiler IRs
- Support for other domains (affine relations, etc.)

But the core algorithm is all here. Go build something!
