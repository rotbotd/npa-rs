# NPA-TP Cheat Sheet

Quick reference for concepts used in npa-rs. Derived from PLAN.md and source.

## Semiring

A set S with two operations:
- **combine** (⊕): like addition, has identity **zero**
- **extend** (⊗): like multiplication, has identity **one**

Laws:
- ⊕ is associative, commutative, has identity 0
- ⊗ is associative, has identity 1, distributes over ⊕
- 0 ⊗ a = a ⊗ 0 = 0

**Kleene star**: a* = 1 ⊕ a ⊕ aa ⊕ aaa ⊕ ... (infinite sum)

## Common Semirings

| Domain | ⊕ | ⊗ | 0 | 1 | a* |
|--------|---|---|---|---|-----|
| Boolean | or | and | false | true | true |
| Shortest path | min | + | ∞ | 0 | 0 |
| Reachability | or | and | false | true | true |
| Sets of paths | ∪ | concat | ∅ | {ε} | Kleene star |
| N×N Bool matrices | pointwise or | matrix mult | zeros | identity | transitive closure |

## The Dataflow Problem

Given a CFG, compute for each node the "summary" of all paths from entry.

For reaching definitions: which definitions might reach this point?
For available expressions: which expressions are definitely computed on all paths?

**Key insight**: This is a fixed-point problem. We want the least solution to:
```
X[entry] = initial
X[n] = ⊕ { f_e(X[m]) | edge e: m → n }
```

## Kleene Iteration

Naive approach: keep iterating until nothing changes.
```
X := 0
repeat:
  X_new := f(X)
until X_new = X
```

**Problem**: Slow for nested loops. Each nesting level can multiply iterations.

## Path Expressions

Instead of iterating, write down a **regular expression** over edges that represents all paths.

Loop `1 → 2 → 3 → 2` with exit `2 → 4`:
- path(2) = a(bc)*
- path(3) = a(bc)*b  
- path(4) = a(bc)*d

Then **evaluate** the regex in your semiring.

## Tarjan's Algorithm

Computes path expressions in nearly-linear time using dominators.

**Dominator**: D dominates N if every path from entry to N goes through D.

Key decomposition:
```
path(N) = path(idom(N)) ⊗ dpath(N)
```

Where dpath(N) = paths from idom(N) to N through nodes dominated by idom(N).

## The Interprocedural Problem

Functions call other functions. Now we have **context** — what happens before and after a call matters.

Equation form:
```
X_j = c_j ⊕ Σ(a ⊗ X_i ⊗ b)
```

The `a ⊗ X ⊗ b` sandwiches are the problem. This is a **Linear Context-Free Language** (LCFL), not regular.

## Newton's Method

From calculus: find roots by linearizing and iterating.

For semirings:
```
ν^(i+1) = f(ν^(i)) ⊕ Df|ν^(i)(Y)
```

Each round solves a **linearized** problem. The linearized version is easier because derivatives are linear.

**Key property**: Converges in ~3-4 rounds typically, vs potentially infinite for Kleene.

## Differentiation Rules

For expressions over semiring:
- D[c] = 0 (constants have zero derivative)
- D[X_j] = Y_j if j = target, else 0
- D[a ⊕ b] = D[a] ⊕ D[b]
- D[a ⊗ b] = D[a] ⊗ b ⊕ a ⊗ D[b] (product rule, but careful about order!)
- D[a*] = a* ⊗ D[a] ⊗ a* (chain rule for star)

The derivative of a* is still **linear in Y** — the a* terms are constants.

## The Non-Commutativity Problem

In normal calculus: d/dx[cx + xd] = c + d

In dataflow, ⊗ doesn't commute! `c ⊗ X + X ⊗ d` stays as two separate terms.

This means linearized equations are still LCFL, not regular. Can't use Tarjan directly.

## Tensor Products (The NPA-TP Trick)

**Admissible semiring** has:
1. **Transpose**: (a ⊗ b)ᵗ = bᵗ ⊗ aᵗ
2. **Tensor product**: S × S → ST (Kronecker product for matrices)
3. **Detensor-transpose**: ST → S, with (t,·)(aᵗ ⊗ b) = a ⊗ b

**The magic**: Transform LCFL to left-linear over tensor semiring, solve with Tarjan, detensor back.

Transform `a ⊗ Y ⊗ b` into `Z ⊗_T (aᵗ ⊗ b)` — now Z is only on the LEFT.

## τ_Reg Transformation

Given LCFL: `Y_j = c_j ⊕ Σ(a ⊗ Y_i ⊗ b)`

Transform to left-linear over ST:
```
Z_j = (1ᵗ ⊗ c_j) ⊕_T Σ(Z_i ⊗_T (aᵗ ⊗ b))
```

Now it's regular! Apply Tarjan in ST, then detensor-transpose to get back to S.

## Coupling Operation

C(a, b) = aᵗ ⊗ b

Key property:
```
C(a₁,b₁) ⊗_T C(a₂,b₂) = C(a₂ ⊗ a₁, b₁ ⊗ b₂)
```

The order reversal in first component + tensor = correct composition!

## For Boolean Matrices (Predicate Abstraction)

- N = 2^|predicates| states
- Elements: N×N Boolean matrices (relations between states)
- ⊕ = union (matrix OR)
- ⊗ = relational composition (matrix multiply with OR/AND)
- Transpose = matrix transpose
- Tensor = Kronecker product (N² × N² matrix)
- Detensor-transpose = existential projection with diagonal constraint

## Algorithm Overview (NPA-TP)

**Preprocessing (once):**
1. Tarjan on each procedure CFG → equations with regular RHS
2. Compute coefficient structure for τ_Reg transformation
3. Build dependence graph between equations
4. Tarjan on dependence graph → parameterized path expressions

**Newton loop:**
```
ν := f(0)
repeat:
  for each (k,j): T_kj := Coeff_k(τReg(D_Xk[Rhs_j]|ν))
  for each j: Z_j := PathExpr_j[substitute T values]
  for each j: ν_j := detensor_transpose(Z_j)
until converged
```

## Implementation Notes

- **Subexpression sharing**: Identical subexpressions should be shared (we use Expr enum with Arc)
- **Memoization**: Cache path expression evaluations (same structure, different leaf values)
- **Sparse matrices**: N×N blows up for large predicate sets, use sparse representation
- **Linear fast path**: For left-linear equations (no sandwiches), skip tensor machinery

## Complexity

- Tarjan: O(m α(m,n)) — essentially linear
- Newton rounds: ~3-4 typically
- Per-round: dominated by matrix operations in semiring

## Why This Beats Kleene

Kleene iteration: each nested loop multiplies iteration count
Newton: fixed number of rounds regardless of nesting depth

The trick: linearization captures the "shape" of recursion in the derivative structure, then we solve the linearized system exactly with path expressions.
