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

