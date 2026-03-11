# Newtonian Program Analysis in Rust

## Background

### The Problem

Interprocedural dataflow analysis requires solving polynomial fixed-point equations over semirings. Kleene iteration converges slowly — logarithmically for probabilistic programs, never terminates for infinite domains.

### Newton's Method for Semirings

Esparza, Kiefer, and Luttenberger (2010) generalized Newton's method to ω-continuous semirings:

```
ν(i+1) = f(ν(i)) ⊕ Df|ν(i)(Y)       # Newton (linear correction)
```

Each round solves a linearized problem. Converges faster than Kleene.

### The Non-Commutativity Problem

In numerical analysis, `c*X + X*d` becomes `(c+d)*X` because multiplication commutes → regular language → Tarjan's path expressions.

In dataflow, "multiplication" is function composition, which does NOT commute. `c*X + X*d` stays as-is → Linear Context-Free Language (LCFL). Harder to solve.

### Tensor Products (Reps et al., POPL 2016)

NPA-TP "regularizes" LCFL problems for admissible semirings:

1. Define tensor product semiring ST (Kronecker product for matrices)
2. Transform LCFL system L into left-linear system LT over ST
3. Solve LT using Tarjan (fast!)
4. Read out via detensor-transpose: `(t,·)(aᵗ ⊗ b) = a ⊗ b`

Key: detensor-transpose distributes over combine — no cross-terms!

---

## Core Concepts

### Admissible Semiring (Definition 4.1)

A semiring S is **admissible** if it has:

1. **Transpose** `·ᵗ : S → S`
   - `(a ⊗ b)ᵗ = bᵗ ⊗ aᵗ` (reverses order)
   - `(a ⊕ b)ᵗ = aᵗ ⊕ bᵗ` (distributes)
   - `(aᵗ)ᵗ = a` (involution)

2. **Tensor product** `⊗ : S × S → ST`
   - Distributes over ⊕ on both sides
   - Mixed product: `(a₁ ⊗ b₁) ⊗T (a₂ ⊗ b₂) = (a₁ ⊗ a₂) ⊗ (b₁ ⊗ b₂)`

3. **Detensor-transpose** `(t,·) : ST → S`
   - `(t,·)(a ⊗ b) = a ⊗ b`
   - `(t,·)(p₁ ⊕T p₂) = (t,·)(p₁) ⊕ (t,·)(p₂)` — DISTRIBUTES!

### Why Cross-Terms Are Avoided

Naive pairing fails because R doesn't distribute over ⊕p:

```
R((a₁,b₁) ⊕p (a₂,b₂)) = (a₁⊕a₂) ⊗ (b₁⊕b₂) = a₁b₁ ⊕ a₂b₁ ⊕ a₁b₂ ⊕ a₂b₂  // CROSS TERMS
```

With tensor products, (t,·) DOES distribute:

```
(t,·)((a₁ᵗ⊗b₁) ⊕T (a₂ᵗ⊗b₂)) = a₁b₁ ⊕ a₂b₂  // NO CROSS TERMS
```

The diagonal constraint in detensor-transpose prevents mixing.

### τ_Reg Transformation (Definition 4.2)

Given LCFL system:

```
Yⱼ = cⱼ ⊕ Σ(aᵢⱼₖ ⊗ Yᵢ ⊗ bᵢⱼₖ)
```

Transform to left-linear over ST:

```
Zⱼ = (1ᵗ ⊗ cⱼ) ⊕T Σ(Zᵢ ⊗T (aᵢⱼₖᵗ ⊗ bᵢⱼₖ))
```

Each `a ⊗ Y ⊗ b` becomes `Z ⊗T (aᵗ ⊗ b)` — left-linear!

### Differential of Kleene-Star (Theorem 6.3)

For loops:

```
D_Xj[g*]|ν(Y) = g(ν)* ⊗ D_Xj[g]|ν(Y) ⊗ g(ν)*
```

Still LINEAR in Y! The g(ν)\* are constants. So we can handle loops directly.

---

## Algorithm 7.1 (NPA-TP)

**Input**: Interprocedural dataflow problem over admissible S, n procedures

**Preprocessing (done ONCE):**

1. **Tarjan on each procedure CFG** → system E with regular RHS (handles loops)

2. **Build coefficient structure**: For each equation Xⱼ, compute `Coeffₖ(τReg(Eⱼ))` = tensor coefficient of Zₖ

3. **Create alphabet symbols** `<k,j>` for the linearized system

4. **Tarjan on linearized system** → parameterized regular expressions `Regⱼ` using alphabet `{<k,j>}`

**Newton Loop:**

```
ν := (f₁(0), ..., fₙ(0))
repeat:
  for k,j = 1 to n:
    Tₖⱼ := Coeffₖ(τReg(Eⱼ))(ν)    // evaluate at current ν
  for j = 1 to n:
    Zⱼ := Regⱼ[<k,j> ← Tₖⱼ]       // substitute and evaluate in ST
    νⱼ := (t,·)(Zⱼ)                // detensor-transpose
until ν converges
```

**Key insight**: The structure of Regⱼ is FIXED. Each round just substitutes new leaf values.

---

## Predicate Abstraction Domain (Section 5)

For Boolean programs with predicate set P:

- N = 2^|P| (number of predicate assignments)
- Semiring elements: N×N Boolean matrices (relations)
- ⊕ = union (matrix OR)
- ⊗ = relational composition (matrix multiply with OR/AND)
- Transpose = matrix transpose
- Tensor = Kronecker product (N² × N² matrix)
- Detensor-transpose = existential projection with diagonal constraint (equation 46):
  ```
  (t,·)(T(A', B, A, B')) = ∃A', B : T(A', B, A, B') ∧ A' = B
  ```

---

## Implementation Plan

### Rust Module Structure

```
src/
├── lib.rs
├── semiring.rs          # Semiring, Admissible traits
├── boolean_matrix.rs    # Predicate abstraction domain
├── tensor.rs            # Kronecker product, detensor-transpose
├── expr.rs              # Regular expressions with variables
├── tarjan.rs            # Path-expression algorithm
├── differentiate.rs     # Differential computation
├── regularize.rs        # τ_Reg transformation
├── npa.rs               # Main NPA-TP loop
└── cfg.rs               # Control flow graph representation
```

### Core Traits

```rust
trait Semiring: Clone + Eq {
    fn zero() -> Self;
    fn one() -> Self;
    fn combine(&self, other: &Self) -> Self;  // ⊕
    fn extend(&self, other: &Self) -> Self;   // ⊗
    fn star(&self) -> Self;                   // Kleene star
}

trait Admissible: Semiring {
    type Tensor: Semiring;

    fn transpose(&self) -> Self;
    fn tensor(&self, other: &Self) -> Self::Tensor;
    fn detensor_transpose(t: &Self::Tensor) -> Self;
}
```

### Implementation Phases

1. **Core Types**
   - [ ] Semiring trait
   - [ ] Boolean matrix type (dense first)
   - [ ] Tensor product type (N² × N²)

2. **Operations**
   - [ ] Matrix combine, extend, star
   - [ ] Kronecker product
   - [ ] Detensor-transpose

3. **Expressions**
   - [ ] Regular expression AST
   - [ ] Generalized regular expression (Defn 4.5)
   - [ ] Expression evaluation

4. **Algorithms**
   - [ ] Tarjan's path-expression algorithm
   - [ ] Differentiation (including Kleene-star)
   - [ ] τ_Reg transformation

5. **NPA-TP Loop**
   - [ ] Coefficient extraction
   - [ ] Parameterized regular expressions
   - [ ] Newton iteration with convergence check

6. **Testing**
   - [ ] Small hand-crafted Boolean programs
   - [ ] Comparison with Kleene iteration
   - [ ] Benchmarks

---

## Open Questions

1. **Sparse representations**: N×N Boolean matrices blow up (N = 2^|predicates|). BDDs? Sparse matrices?

2. **Incremental updates**: When program changes, reuse Tarjan structure?

3. **Beyond predicate abstraction**: What other semirings are admissible? (affine relations mentioned)

4. **Widening**: For infinite-height lattices, how does widening interact with Newton?

---

## Algorithm 7.1 (Full Detail from Paper)

**Input**: Interprocedural dataflow problem over admissible S, procedures X̃

1. **Tarjan on each procedure CFG** → system E of equations with regular RHS
   - Each variable = one procedure
   - RHS = regular expression over variables and constants

2. **For each equation Xⱼ = Rhsⱼ(X̃)**, create left-linear equation:

   ```
   Zⱼ = τReg(D Rhsⱼ|ν̃(Ỹ))
   ```

3. **Create dependence graph G**:
   - Edge Zₖ → Zⱼ labeled ⟨k,j⟩ if Zⱼ's equation contains Zₖ on RHS
   - Dummy vertex Λ with edges Λ → Zⱼ labeled ⟨0,j⟩ for each Zⱼ

4. **Tarjan on G** (with entry Λ) → regular expression Rᵢ for each Zᵢ
   - Alphabet: {⟨k,j⟩ | 0 ≤ k ≤ n, 1 ≤ j ≤ n}

5. **Create map m**: variable Zᵢ maps to:

   ```
   Rᵢ[⟨0,j⟩ ← (1ᵗ ⊗ Rhsⱼ(ν̃))]
     [⟨1,j⟩ ← Coeff₁(τReg(D_X₁ Rhsⱼ|ν̃(Ỹ)))]
     ...
     [⟨n,j⟩ ← Coeffₙ(τReg(D_Xₙ Rhsⱼ|ν̃(Ỹ)))]
   ```

6. **Initialize**: i ← 0; μ̃ ← f̃(0̃)

7. **Repeat**:
   - ν̃⁽ⁱ⁾ = μ̃
   - μ̃ = ⟨(t,·)([[m(Zⱼ)]]\_T ν̃⁽ⁱ⁾) | Zⱼ ∈ Z̃⟩
   - i ← i + 1

   **until** ν̃⁽ⁱ⁻¹⁾ = μ̃

8. **Return** μ̃

## Local Variables (Section 8)

For programs with local variables, use merge functions:

**Merge function** M : D × D → D satisfies:

- 0-strictness: M(a, 0) = 0, M(0, b) = 0
- Distributivity: M distributes over ⊕ in both positions
- Path extension: M(a ⊗ b, c) = a ⊗ M(b, c)

**Project operator**: Project(R) = M(1, R)

For relational weight domain with globals G and locals L:

```
Project(R(G, L, G', L')) = (∃L, L' : R) ∧ (L = L')
```

For tensor-product semiring:

```
ProjectT(T(G'₁, L'₁, G₂, L₂, G₁, L₁, G'₂, L'₂)) =
  (∃L'₁, L₂, L₁, L'₂ : T ∧ (L'₁ = L₂)) ∧ (L₁ = L'₁) ∧ (L₂ = L'₂)
```

**Differential of Project**:

```
D_Xj[Project(g)]|ν(Y) = Project(D_Xj[g]|ν(Y))
```

When using NPA-TP with local variables, insert ProjectT after step 2:

```
Rᵢ[⟨0,j⟩ ← ProjectT(1ᵗ ⊗ Rhsⱼ(ν̃))]
  [⟨k,j⟩ ← ProjectT(Coeffₖ(τReg(D_Xₖ Rhsⱼ|ν̃(Ỹ))))]
```

## Experimental Results (Section 9)

From POPL16 experiments on 584 Boolean programs:

| Solver | Completed | Timeouts | Spaceouts | Avg Newton Rounds |
| ------ | --------- | -------- | --------- | ----------------- |
| EWPDS  | 495       | 16       | 73        | N/A               |
| FWPDS  | 483       | 32       | 69        | N/A               |
| NPA    | 290       | 142      | 152       | 3.38              |
| NPA-TP | 386       | 16       | 182       | 3.67              |

Key findings:

- NPA-TP faster than NPA (geometric means: 1.62 → 4.61)
- EWPDS still faster than NPA-TP (geometric means: 5.20 → 6.84)
- Newton converges in ~3-4 rounds on average
- For this test suite, chaotic iteration (EWPDS) won overall

**But**: NPA-TP may win for domains where Kleene-star is expensive (non-idempotent semirings, probabilistic programs).

## Tarjan's Path-Expression Algorithm

Two uses in NPA-TP:

1. **Step 1**: Applied to each procedure CFG to create equations with regular RHS (handles loops via Kleene-star)

2. **Step 4**: Applied to dependence graph G to create parameterized regular expressions

Tarjan's algorithm: O(m α(m,n)) where α is inverse Ackermann. Can degenerate on non-reducible graphs.

## Implementation Notes

- **Subexpression sharing**: Identical subexpressions of regular expressions should be shared
- **Function caching**: Use memoization during regular-expression evaluation (step 7b)
- **OBDD variable ordering**: For predicate abstraction, OBDD size is sensitive to variable ordering

## Tensor-Product Principle (from Related Work)

> Tensor products—plus an appropriate detensor operation—allow computations to be rearranged in certain ways; they can be used to delay performing every multiplication in a sequence of multiplications, which is useful if either:
> (a) a value that is only obtainable at a later time needs to be placed in the middle of the sequence, or
> (b) a subsequence of values in the middle of the sequence needs to be adjusted in certain ways before contributing to the overall product.

## Generalized Regular Expressions (Definition 4.5)

The regular expressions in NPA-TP involve mixed operators:

```
expT ::= aT ∈ ST
       | expᵗ ⊗ exp        // coupling: creates ST from S
       | expT ⊕T expT
       | expT ⊗T expT
       | expT *T

exp  ::= a ∈ S
       | νi               // variable symbols (current round values)
       | exp ⊕ exp
       | exp ⊗ exp
       | exp *
```

Evaluation with respect to value vector ν̃:

- Constants evaluate to themselves
- νi evaluates to (ν̃)i (lookup)
- Operators interpreted in S or ST as appropriate
- Coupling `aᵗ ⊗ b` creates ST value from S values

## Local Variables (Section 8)

For programs with local variables (globals G, locals L):

**Merge function** M : D × D → D satisfies:

- 0-strictness: M(a, 0) = 0, M(0, b) = 0
- Distributivity over ⊕ in both positions
- Path extension: M(a ⊗ b, c) = a ⊗ M(b, c)

**Project operator**: Project(R) = M(1, R)

For relational domain with G and L:

```
Project(R(G, L, G', L')) = (∃L, L' : R) ∧ (L = L')
```

For tensor semiring:

```
ProjectT(T(G'₁, L'₁, G₂, L₂, G₁, L₁, G'₂, L'₂)) =
  (∃L'₁, L₂, L₁, L'₂ : T ∧ (L'₁ = L₂)) ∧ (L₁ = L'₂) ∧ (L₂ = L'₂)
```

**Key properties of ProjectT**:

- ProjectT(ProjectT(a) ⊗T b) = ProjectT(a) ⊗T ProjectT(b)
- ProjectT(a ⊕T b) = ProjectT(a) ⊕T ProjectT(b)
- ProjectT(ProjectT(a)) = ProjectT(a)

**Differential of Project**:

```
D_Xj[Project(g)]|ν(Y) = Project(D_Xj[g]|ν(Y))
```

When using NPA-TP with locals, insert ProjectT in step 5:

```
Rᵢ[⟨0,j⟩ ← ProjectT(1ᵗ ⊗ Rhsⱼ(ν̃))]
  [⟨k,j⟩ ← ProjectT(Coeffₖ(τReg(D_Xₖ Rhsⱼ|ν̃(Ỹ))))]
```

## Experimental Results (Section 9)

From POPL16 on 584 Boolean programs:

| Solver | Completed | Timeouts | Spaceouts | Avg Rounds |
| ------ | --------- | -------- | --------- | ---------- |
| EWPDS  | 495       | 16       | 73        | N/A        |
| FWPDS  | 483       | 32       | 69        | N/A        |
| NPA    | 290       | 142      | 152       | 3.38       |
| NPA-TP | 386       | 16       | 182       | 3.67       |

Findings:

- NPA-TP faster than NPA (geomean: 1.62→4.61)
- EWPDS still faster than NPA-TP (geomean: 5.20→6.84)
- Newton converges in ~3-4 rounds
- For Boolean programs, chaotic iteration (EWPDS) won overall

**But**: NPA-TP may win for non-idempotent semirings (probabilistic programs).

## Implementation Notes (from paper)

1. **Subexpression sharing**: Identical subexpressions should be shared
2. **Function caching**: Use memoization during regular-expression evaluation
3. **OBDD variable ordering**: For predicate abstraction, OBDD size is sensitive to ordering
4. **Tarjan degeneracy**: Can degenerate on non-reducible graphs (didn't happen in experiments)

## Tensor-Product Principle (from Related Work)

> Tensor products—plus an appropriate detensor operation—allow computations to be rearranged in certain ways; they can be used to delay performing every multiplication in a sequence of multiplications, which is useful if either:
> (a) a value that is only obtainable at a later time needs to be placed in the middle of the sequence, or
> (b) a subsequence of values in the middle of the sequence needs to be adjusted in certain ways before contributing to the overall product.

## References

- Esparza, Kiefer, Luttenberger. "Newtonian Program Analysis." JACM 2010.
- Reps, Turetsky, Prabhu. "Newtonian Program Analysis via Tensor Product." POPL 2016.
- Tarjan. "Fast algorithms for solving path problems." JACM 1981.
- Lal, Reps. "Improving pushdown system model checking." CAV 2006. (FWPDS)
- Lal, Reps, Balakrishnan. "Extended weighted pushdown systems." CAV 2005. (EWPDS)

---

## Algorithm 3.4 (Naive Pairing - What Fails)

The naive approach that DOESN'T work, but motivates tensor products:

1. Convert LCFL system L into left-linear system LReg using paired values
2. Solve LReg
3. Apply readout R to get solution to L
   **Why it fails**: R doesn't distribute over ⊕p, causing cross-terms.
   The paper shows 8 summands arise from R(AB), but only 4 meet the matching condition. The other 4 are cross-terms.

---

## Example from Paper (Equation 3)

Original equation: `X1 = a(cX1d)*b`
After linearization:

```
Y1 = a(cν1d)*b ⊕ a(cν1d)* Y1 (cν1d)*b
```

After τ_Reg:

```
Z1 = (1ᵗ ⊗ a(cν1d)*b) ⊕T Z1 ⊗T ((a(cν1d)*)ᵗ ⊗ (cν1d)*b)
```

Regular expression for Z1: `A ⊗T B*T` where:

- A = (1ᵗ ⊗ (d ⊕ bν2ν2c))
- B = (bᵗ ⊗ (ν2c)) ⊕T ((bν2)ᵗ ⊗ c)

---

## Full Detensor-Transpose Derivation (Section 4.4)

The paper walks through why (t,·)(AB) produces exactly 4 correct terms:

```
(t,·)(AB) = (t,·)((1ᵗ ⊗ (d⊕bν2ν2c)) ⊗T ((bᵗ ⊗ ν2c) ⊕T ((bν2)ᵗ ⊗ c)))
```

Using distributivity of ⊗T over ⊕T:

```
= (t,·)((1ᵗ ⊗ (d⊕bν2ν2c)) ⊗T (bᵗ ⊗ ν2c))
  ⊕ (t,·)((1ᵗ ⊗ (d⊕bν2ν2c)) ⊗T ((bν2)ᵗ ⊗ c))
```

Using mixed-product property (a1⊗b1) ⊗T (a2⊗b2) = (a1⊗a2) ⊗ (b1⊗b2):

```
= (t,·)((1ᵗ ⊗ bᵗ) ⊗ ((d⊕bν2ν2c) ⊗ ν2c))
  ⊕ (t,·)((1ᵗ ⊗ (bν2)ᵗ) ⊗ ((d⊕bν2ν2c) ⊗ c))
```

Simplifying and applying (t,·):

```
= bdν2c ⊕ bbν2ν2cν2c ⊕ bν2dc ⊕ bν2bν2ν2cc
```

## Exactly the 4 correct terms! The position of "⊗" marks where the recursive substitution happens.

## Key Implementation Details

### Subexpression Sharing

Identical subexpressions of regular expressions should be shared. This is critical for efficiency.

### Function Caching

Use memoization during regular-expression evaluation (step 7b). Same structure evaluated with different ν values.

### OBDD Variable Ordering

For predicate abstraction with OBDDs, size is very sensitive to variable ordering. This affects performance significantly.

### Tarjan Degeneracy

## Can degenerate on non-reducible CFGs, but didn't happen in experiments.

## Coupling Operation (Definition)

The coupling of a and b is:

```
C(a, b) = aᵗ ⊗ b
```

Key property (equation 35):

```
C(a1,b1) ⊗T C(a2,b2) = (a1ᵗ ⊗ b1) ⊗T (a2ᵗ ⊗ b2)
                      = (a1ᵗ ⊗ a2ᵗ) ⊗ (b1 ⊗ b2)
                      = (a2 ⊗ a1)ᵗ ⊗ (b1 ⊗ b2)
                      = C(a2 ⊗ a1, b1 ⊗ b2)
```

And detensor-transpose recovers:

```
(t,·)(C(a2⊗a1, b1⊗b2)) = a2 ⊗ a1 ⊗ b1 ⊗ b2
```

## The order reversal in the first component + tensor product = correct pairing!

## Why NPA-TP May Still Win

From the paper's conclusion:

> NPA-TP is still slower than EWPDS... However, for domains where Kleene-star is expensive (non-idempotent semirings, probabilistic programs), NPA-TP may be the algorithm of choice.
> The experiments used Boolean programs (idempotent semiring). For:

- Probabilistic programs
- Affine relations
- Other non-idempotent domains
  NPA-TP's faster convergence (3-4 rounds) may outweigh the per-round overhead.

---

## Tarjan's Path-Expression Algorithm

From: https://rolph-recto.github.io/blog/introduction-to-tarjans-path-expressions-algorithm.html

### The Problem

Computing path expressions = regular expressions over edges representing all paths from entry to each node.

Kleene's algorithm is O(n³) and computes ALL pairs — overkill when we only need entry → all nodes.

### Key Insight: Decompose with Dominators

**Dominator**: Node D dominates N if every path from entry to N goes through D.

**Immediate dominator**: idom(N) = the unique closest strict dominator of N.

**Key lemma**: For any path to N, after the last occurrence of some dominator D, all remaining nodes are properly dominated by D. Paths "go down" the dominator tree.

### Decomposition

```
path(N) = path(idom(N)) × dpath(N)
```

Where `dpath(N)` = paths from idom(N) to N that only visit nodes dominated by idom(N).

Recursively:

```
path(N) = path(D₁) × dpath(D₂) × ... × dpath(Dₖ)
```

Where D₁...Dₖ is the path from entry to N in the dominator tree.

### Computing dpath

Post-order traversal of dominator tree. When processing node N:

1. **Tree edges**: `tree(C)` = incoming edges to child C from N (immediate dominator)

   ```
   dpathTree(C) = e₁ + e₂ + ... + eₙ   (sum of tree edges)
   ```

2. **Non-tree edges**: `nontree(C)` = other incoming edges to C
   - Source of non-tree edge must be dominated by idom(target) (Lemma 2)
   - So source is dominated by some sibling C' of C in dominator tree

3. **Derived graph**: For each non-tree edge e to C, find sibling C' that dominates source(e), add edge C' → C in derived graph

4. **Image of derived edge**:

   ```
   image(e') = dpath(D₂) × ... × dpath(Dₖ) × dpath(source(e)) × e
   ```

   Where D₁...Dₖ,source(e) is path from C' to source(e) in domtree

5. **Compute dpath**:
   ```
   dpath(C) = Σ_{C' ∈ siblings(C)} dpathTree(C') × img(derivedPath(C', C))
   ```

### Computing path(entry)

Entry has no immediate dominator. For each backedge e to entry:

```
Pₑ = dpath(D₂) × ... × dpath(Dₖ) × e
```

Where D₁...Dₖ,source(e) is path from entry to source(e) in domtree.

Then:

```
path(entry) = (Σₑ Pₑ)*
```

### Example

CFG with loop: 1 → 2 → 3 → 2, and 2 → 4

Dominator tree: 1 dominates 2, 2 dominates {3, 4}

Derived graph: edge c' from 2 to itself (back edge 3→2)

Computing:

- dpath(3) = b
- dpath(4) = d
- image(c') = dpath(3) × c = bc
- derivedPath(2,2) = (c')\*
- dpath(2) = dpathTree(2) × img((c')_) = a × (bc)_ = a(bc)\*

Results:

- path(1) = 1
- path(2) = a(bc)\*
- path(3) = a(bc)\*b
- path(4) = a(bc)\*d

### Complexity

O(m α(m,n)) where α is inverse Ackermann — essentially linear.

The key: decomposition via dominators breaks the problem into smaller pieces that can be solved independently.

### Implementation Notes

1. Compute dominator tree first (Lengauer-Tarjan algorithm)
2. Post-order traverse domtree to compute dpath
3. Reverse post-order to compute path (composing with dpath)
4. For derived graph path expressions between siblings, use simpler algorithm (Kleene-like "eliminate")

Python implementation: https://gist.github.com/rolph-recto/0e1ffb042f690fd09387920d24677ba4

### Key Functions for Implementation

```python
def tarjan(cfg, domtree):
    # derivedEdges = edges of derived graph
    # edgeMap = maps CFG edges to derived graph edges
    derivedEdges, edgeMap = computeDerivedGraph(cfg, domtree)

    path = dict()
    dpath = dict()
    dpathTree = dict()

    # compute dpath in POST-ORDER traversal of domtree
    for node in postorder(domtree):
        children = domtree.children(node)
        imageMap = dict()

        # compute dpathTree and image of derived edges
        for child in children:
            dpathTree[child] = sum(tree(child))  # sum of tree edges
            for edge in nontree(child):
                siblingDom = edgeMap[edge].src
                imageMap[edgeMap[edge]] = edgePaths(edge, siblingDom, dpath, domtree)

        # compute path expressions in derived graph between siblings
        derivedPath = computePathExpressions(children, derivedEdges)

        # compute dpath for each child
        for child in children:
            dpath[child] = sum(
                dpathTree[sibling] * pathImg(imageMap, derivedPath[(sibling, child)])
                for sibling in children
            )

    # compute path in REVERSE POST-ORDER traversal
    for node in reverse_postorder(domtree):
        if node == cfg.entry:
            # special case: entry has no idom
            path[entry] = star(sum(edgePaths(e, entry, dpath, domtree) for e in nontree(entry)))
        else:
            path[node] = path[idom(node)] * dpath[node]

    return path
```

### edgePaths Function

For an edge `e` from `source(e)` to `target(e)`, and a starting node `start`:

```
edgePaths(e, start, dpath, domtree) = dpath(D₂) × ... × dpath(Dₖ) × dpath(source(e)) × e
```

Where `D₁, D₂, ..., Dₖ, source(e)` is the path from `start` to `source(e)` in the dominator tree.

---

## Algorithm 3.4 (Naive Pairing - What Fails)

The naive approach that DOESN'T work, but motivates tensor products:

1. Convert LCFL system L into left-linear system LReg using paired values
2. Solve LReg via Tarjan
3. Apply readout R(a,b) = a ⊗ b to get solution to L

**Paired semiring operations:**

- `(a₁, b₁) ⊗p (a₂, b₂) = (a₂ ⊗ a₁, b₁ ⊗ b₂)` — note order reversal!
- `(a₁, b₁) ⊕p (a₂, b₂) = (a₁ ⊕ a₂, b₁ ⊕ b₂)`

**Why it fails**: R doesn't distribute over ⊕p, causing cross-terms.

The paper shows 8 summands arise from R(AB), but only 4 meet the matching condition. The other 4 are cross-terms from `(a₁⊕a₂) ⊗ (b₁⊕b₂)` expanding.

---

## Coupling Operation (Definition)

The coupling of a and b is:

```
C(a, b) = aᵗ ⊗ b
```

Key property (equation 35):

```
C(a₁,b₁) ⊗T C(a₂,b₂) = (a₁ᵗ ⊗ b₁) ⊗T (a₂ᵗ ⊗ b₂)
                      = (a₁ᵗ ⊗ a₂ᵗ) ⊗ (b₁ ⊗ b₂)
                      = (a₂ ⊗ a₁)ᵗ ⊗ (b₁ ⊗ b₂)
                      = C(a₂ ⊗ a₁, b₁ ⊗ b₂)
```

And detensor-transpose recovers:

```
(t,·)(C(a₂⊗a₁, b₁⊗b₂)) = a₂ ⊗ a₁ ⊗ b₁ ⊗ b₂
```

The order reversal in the first component + tensor product = correct pairing!

---

## Full Detensor-Transpose Derivation (Section 4.4)

The paper walks through why (t,·)(AB) produces exactly 4 correct terms.

Key steps:

1. Distribute ⊗T over ⊕T (equation 40→41)
2. Apply mixed-product property: `(a₁⊗b₁) ⊗T (a₂⊗b₂) = (a₁⊗a₂) ⊗ (b₁⊗b₂)`
3. Simplify and apply (t,·)
4. Result: exactly the 4 terms that satisfy the matching condition

The position of "⊗" in `bᵗ ⊗ (ν₂⊗c)` marks where the recursive substitution happens.

---

## Tensor-Product Principle (from Related Work)

> Tensor products—plus an appropriate detensor operation—allow computations to be rearranged in certain ways; they can be used to delay performing every multiplication in a sequence of multiplications, which is useful if either:
> (a) a value that is only obtainable at a later time needs to be placed in the middle of the sequence, or
> (b) a subsequence of values in the middle of the sequence needs to be adjusted in certain ways before contributing to the overall product.

---

## Why NPA-TP May Still Win

From the paper's conclusion:

> NPA-TP is still slower than EWPDS... However, for domains where Kleene-star is expensive (non-idempotent semirings, probabilistic programs), NPA-TP may be the algorithm of choice.

The experiments used Boolean programs (idempotent semiring). For:

- Probabilistic programs
- Affine relations
- Other non-idempotent domains
