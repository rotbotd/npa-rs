# Lecture 7: Going Interprocedural with Tensor Products

Last lecture we saw that Newton's method linearizes the problem but the linearized system is still LCFL. This lecture shows how tensor products rescue us.

## The Setup

After Newton linearization, we have equations like:

```
Y_j = c_j ⊕ Σ(a ⊗ Y_i ⊗ b)
```

The `a ⊗ Y ⊗ b` sandwiches are the problem. We need to eliminate them.

## Naive Pairing (What Doesn't Work)

First, let's understand what we're trying to do and why the obvious approach fails.

Idea: Instead of tracking just the left context or just the right context, track **both** as a pair (a, b). Then `a ⊗ Y ⊗ b` becomes `Y ⊗ (a, b)` — the unknown is on the left!

Define a paired semiring:
- `(a₁, b₁) ⊗p (a₂, b₂) = (a₂ ⊗ a₁, b₁ ⊗ b₂)` — note the order reversal!
- `(a₁, b₁) ⊕p (a₂, b₂) = (a₁ ⊕ a₂, b₁ ⊕ b₂)`

Readout: `R(a, b) = a ⊗ b`

**Why it fails**: R doesn't distribute over ⊕p!

```
R((a₁,b₁) ⊕p (a₂,b₂)) = R(a₁⊕a₂, b₁⊕b₂)
                       = (a₁⊕a₂) ⊗ (b₁⊕b₂)
                       = a₁⊗b₁ ⊕ a₁⊗b₂ ⊕ a₂⊗b₁ ⊕ a₂⊗b₂
```

But we wanted:

```
R(a₁,b₁) ⊕ R(a₂,b₂) = a₁⊗b₁ ⊕ a₂⊗b₂
```

The extra terms `a₁⊗b₂` and `a₂⊗b₁` are **cross-terms**. They don't correspond to any actual path — they mix the left context from one path with the right context from another.

## Tensor Products to the Rescue

The tensor product construction keeps pairs apart by living in a larger space where cross-terms can't form.

### The Tensor Product Semiring

For a semiring S, the tensor product semiring ST has:
- Elements: formal sums of pairs `aᵗ ⊗ b` where a, b ∈ S
- `⊕T`: union of formal sums
- `⊗T`: the **mixed product** — `(a₁ᵗ ⊗ b₁) ⊗T (a₂ᵗ ⊗ b₂) = (a₁ᵗ ⊗ a₂ᵗ) ⊗ (b₁ ⊗ b₂)`

For Boolean matrices, this is the **Kronecker product** — if S is N×N matrices, ST is N²×N² matrices.

### Detensor-Transpose

The magic operation: `(t,·) : ST → S`

```
(t,·)(aᵗ ⊗ b) = a ⊗ b
```

**Key property**: (t,·) distributes over ⊕T!

```
(t,·)(p₁ ⊕T p₂) = (t,·)(p₁) ⊕ (t,·)(p₂)
```

No cross-terms! This is because the tensor structure keeps the pairs separate until we explicitly detensor.

### Admissible Semirings

A semiring is **admissible** if it has:

1. **Transpose**: `·ᵗ : S → S`
   - `(a ⊗ b)ᵗ = bᵗ ⊗ aᵗ` (reverses order)
   - `(a ⊕ b)ᵗ = aᵗ ⊕ bᵗ`
   - `(aᵗ)ᵗ = a`

2. **Tensor product**: `⊗ : S × S → ST`
   - Creates the tensor semiring

3. **Detensor-transpose**: `(t,·) : ST → S`
   - Distributes over ⊕T

Boolean matrices are admissible — transpose is matrix transpose, tensor is Kronecker product.

## The τ_Reg Transformation

Given an LCFL system:

```
Y_j = c_j ⊕ Σ(a ⊗ Y_i ⊗ b)
```

Transform to left-linear over ST:

```
Z_j = (1ᵗ ⊗ c_j) ⊕T Σ(Z_i ⊗T (aᵗ ⊗ b))
```

**The key move**: `a ⊗ Y ⊗ b` becomes `Z ⊗T (aᵗ ⊗ b)`

- Before: Y is sandwiched between a and b
- After: Z is on the LEFT, the context `(aᵗ ⊗ b)` is on the right

This is left-linear! We can use Tarjan!

## The Coupling Operation

Define **coupling**: `C(a, b) = aᵗ ⊗ b`

Coupling has a beautiful composition property:

```
C(a₁,b₁) ⊗T C(a₂,b₂) = C(a₂ ⊗ a₁, b₁ ⊗ b₂)
```

Notice:
- Left components compose in **reverse** order: a₂ ⊗ a₁
- Right components compose in **forward** order: b₁ ⊗ b₂

This is exactly what we need for function composition! When you go through function f then function g:
- The pre-context accumulates as g's then f's (reverse)
- The post-context accumulates as f's then g's (forward)

And detensor-transpose recovers the final composition:

```
(t,·)(C(a₂⊗a₁, b₁⊗b₂)) = (a₂⊗a₁) ⊗ (b₁⊗b₂)
```

## Worked Example

Original equation: `X = a ⊗ X ⊗ b ⊕ c`

After τ_Reg: `Z = Z ⊗T (aᵗ ⊗ b) ⊕T (1ᵗ ⊗ c)`

This is left-linear: `Z = Z ⊗T A ⊕T C` where:
- `A = aᵗ ⊗ b` (tensor coefficient)
- `C = 1ᵗ ⊗ c` (constant term)

Solution: `Z = C ⊗T A* = (1ᵗ ⊗ c) ⊗T (aᵗ ⊗ b)*`

Detensor: `X = (t,·)(Z)`

## For Boolean Matrices

If S is N×N Boolean matrices:
- ST is N²×N² Boolean matrices
- Transpose: `Aᵗ[i,j] = A[j,i]`
- Tensor (Kronecker): `(A ⊗ B)[(i,j), (k,l)] = A[i,k] ∧ B[j,l]`
- Detensor-transpose: existential projection with diagonal constraint

The detensor-transpose formula (equation 46 in paper):

```
(t,·)(T)[A, A'] = ∃B, B' : T[(A',B), (A,B')] ∧ (B = B')
```

In words: T is a relation on pairs (pre, post). Detensor gives the relation where the "inner" post (B) equals the "inner" pre (B').

## The Full NPA-TP Algorithm

Putting it all together:

**Preprocessing (once):**
1. Tarjan on each procedure CFG → equations with regular RHS
2. Compute coefficient structure for τ_Reg
3. Tarjan on dependence graph → parameterized path expressions

**Newton loop:**
```
ν := f(0)
repeat:
  // Compute tensor coefficients at current ν
  for each (k,j):
    T_kj := Coeff_k(τReg(D_Xk[Rhs_j]|ν))
  
  // Solve left-linear system over ST
  for each j:
    Z_j := PathExpr_j with T values substituted
  
  // Detensor to get new ν
  for each j:
    ν_j := (t,·)(Z_j)
    
until ν converges
```

**Key insight**: The path expression structure is fixed! Each round just substitutes new leaf values.

## Why This Works

1. **Newton** reduces nonlinear → linear in ~3-4 rounds
2. **τ_Reg** reduces LCFL → left-linear over tensor semiring
3. **Tarjan** solves left-linear in nearly-linear time
4. **Detensor** recovers the solution without cross-terms

Each piece handles one source of complexity. Together: polynomial time for what would otherwise be exponential.

## In Code

From our implementation:

```rust
pub trait Admissible: Semiring {
    type Tensor: Semiring;

    fn transpose(&self) -> Self;
    fn tensor(&self, other: &Self) -> Self::Tensor;
    fn detensor_transpose(t: &Self::Tensor) -> Self;
}

/// Coupling operation: C(a, b) = aᵗ ⊗ b
pub fn coupling<S: Admissible>(a: &S, b: &S) -> S::Tensor {
    a.transpose().tensor(b)
}
```

## Problem Set 7

### Problem 7.1: Cross-Terms

For the naive pairing approach with:
- `(a₁, b₁) = (X, Y)`
- `(a₂, b₂) = (Z, W)`

Compute:
a) `(a₁, b₁) ⊕p (a₂, b₂)`
b) `R((a₁, b₁) ⊕p (a₂, b₂))` where R(a,b) = a ⊗ b
c) `R(a₁, b₁) ⊕ R(a₂, b₂)`
d) What are the cross-terms?

### Problem 7.2: Kronecker Product

For 2×2 Boolean matrices:
```
A = [[1,0],[0,1]]  (identity)
B = [[1,1],[0,0]]  (first row all 1s)
```

Compute the Kronecker product `A ⊗ B` (a 4×4 matrix).

### Problem 7.3: Coupling Composition

Verify the coupling composition property for 2×2 Boolean matrices:

```
C(a₁,b₁) ⊗T C(a₂,b₂) = C(a₂ ⊗ a₁, b₁ ⊗ b₂)
```

Use:
- `a₁ = [[1,0],[0,0]]`, `b₁ = [[0,0],[0,1]]`
- `a₂ = [[1,1],[0,0]]`, `b₂ = [[1,0],[1,0]]`

### Problem 7.4: τ_Reg by Hand

Transform this LCFL equation using τ_Reg:

```
Y = a ⊗ Y ⊗ b ⊕ c ⊗ Y ⊗ d ⊕ e
```

Write the resulting left-linear equation over the tensor semiring.

### Problem 7.5: Detensor-Transpose Semantics

For Boolean matrices representing relations:
- State space: {0, 1}
- A tensor element T is a 4×4 Boolean matrix over pairs {(0,0), (0,1), (1,0), (1,1)}

If T[(0,0), (1,1)] = true and all else false, what is (t,·)(T)?

Interpret this in terms of relations.

### Problem 7.6: The Diagonal Constraint

The detensor-transpose formula has a diagonal constraint `B = B'`.

Why is this necessary? What would go wrong if we just projected without the constraint?

Hint: Think about what happens at a function call — the callee's entry state should match the call site's post-state.

### Problem 7.7: Size Blowup

For predicate abstraction with k predicates:
- Base matrices: 2^k × 2^k
- Tensor matrices: 2^(2k) × 2^(2k)

With k = 10 predicates:
- How big is the base matrix?
- How big is the tensor matrix?

This is why sparse representations matter!

### Problem 7.8: Solving a Small Example

Consider the mutual recursion:
```
X = a ⊗ Y ⊗ b
Y = c ⊗ X ⊗ d ⊕ e
```

Using τ_Reg, write the transformed system over ST.

Then argue (without computing) how many Newton rounds you'd expect before convergence.

---

Next: [Lecture 8: Implementing It All](./08-implementation.md)
