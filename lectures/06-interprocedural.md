# Lecture 6: The Interprocedural Problem

We have path expressions. We have Tarjan's algorithm. Intraprocedural dataflow is solved. But what about function calls?

## Intraprocedural: Solved

For a single procedure with loops, path expressions work beautifully. The loop becomes a Kleene star, and we evaluate it directly.

```
1 → 2 → 3 → 2 (back)
    ↓
    4
```

Path expression for node 4: `a ⊗ (bc)* ⊗ d`

Evaluate in whatever semiring you want. Done.

## The Interprocedural Problem

Now add function calls:

```
proc main:
    1: call foo
    2: x = result
    3: call foo
    4: y = result

proc foo:
    5: if (...) call bar else return 0
    6: return result + 1

proc bar:
    7: return 42
```

Each procedure has a **summary**: what does it do to the dataflow state? The summary of `main` depends on the summary of `foo`. The summary of `foo` depends on `bar`. If they were mutually recursive, they'd depend on each other.

We need to solve:

```
X_main = ... f(X_foo) ...
X_foo  = ... g(X_bar) ...
X_bar  = ... h() ...
```

This is a system of **polynomial equations** over the semiring.

## The Sandwich Problem

Here's where it gets ugly. A procedure call looks like:

```
(stuff before call) ⊗ X_callee ⊗ (stuff after call)
```

The callee's summary is **sandwiched** between the caller's context.

Written as an equation:

```
X = a ⊗ Y ⊗ b ⊕ c
```

The variable Y appears in the **middle** of a product. Not on the left, not on the right — in the middle.

## Why This Breaks Path Expressions

Path expressions work when everything is **left-linear** or **right-linear**:

- Left-linear: `Y = a ⊗ Y ⊕ c` → solve with `Y = a* ⊗ c`
- Right-linear: `Y = Y ⊗ b ⊕ c` → solve with `Y = c ⊗ b*`

But sandwiches are **bilinear**:

```
Y = a ⊗ Y ⊗ b ⊕ c
```

If ⊗ were commutative, we could rewrite: `a ⊗ Y ⊗ b = (a ⊗ b) ⊗ Y`

But ⊗ is NOT commutative for matrices! Relational composition doesn't commute. `(A then B) ≠ (B then A)`.

This is a **Linear Context-Free Language** (LCFL), not a regular language. Tarjan doesn't apply directly.

## What Happens With Kleene Iteration?

You might think: just iterate. We know iteration converges.

It does... eventually. But for nested structures, convergence is SLOW.

Consider mutual recursion:

```
X = a ⊗ Y ⊗ b
Y = c ⊗ X ⊗ d ⊕ e
```

Unfolding:
```
X = a ⊗ (c ⊗ X ⊗ d ⊕ e) ⊗ b
  = a ⊗ c ⊗ X ⊗ d ⊗ b ⊕ a ⊗ e ⊗ b
```

X depends on X, sandwiched. Each iteration adds another layer of context. For k levels of mutual recursion, you might need O(2^k) iterations.

## The Insight: Newton's Method

In calculus, Newton's method finds roots by repeatedly linearizing:

```
f(x) = 0
x_{n+1} = x_n - f(x_n) / f'(x_n)
```

Each step solves a **linear** approximation to the nonlinear problem.

For semirings, Esparza, Kiefer, and Luttenberger (2010) showed:

```
ν^{(i+1)} = f(ν^{(i)}) ⊕ Df|_{ν^{(i)}}(Y)
```

Where `Df|_ν(Y)` is the **formal derivative** of f with respect to the variables, evaluated at ν.

The linearized problem is easier! It converges in ~3-4 rounds typically, regardless of nesting depth.

## Derivatives of Semiring Expressions

For an expression over a semiring, the derivative rules are:

- `D[c] = 0` (constants have zero derivative)
- `D[X_j]` with respect to X_k: `Y_k` if j=k, else `0`
- `D[a ⊕ b] = D[a] ⊕ D[b]`
- `D[a ⊗ b] = D[a] ⊗ b ⊕ a ⊗ D[b]` (product rule!)
- `D[a*] = a* ⊗ D[a] ⊗ a*` (chain rule for star)

Notice: the derivative of a* involves **a* on both sides of D[a]**. But those are **constants** (evaluated at the current ν). So the derivative is still **linear** in the new variable Y.

## The Linearized System Is Still LCFL

Here's the problem. After differentiation:

```
D[a ⊗ X ⊗ b] = D[a] ⊗ X ⊗ b ⊕ a ⊗ D[X] ⊗ b ⊕ a ⊗ X ⊗ D[b]
```

The middle term has `D[X] = Y` sandwiched between `a` and `b`. We STILL have a sandwich!

The linearized system is:

```
Y = ... ⊕ a ⊗ Y ⊗ b ⊕ ...
```

This is still LCFL, not regular. We can't use Tarjan directly on the linearized system.

## Preview: Tensor Products

This is where tensor products come in (Lecture 7). The key insight:

Transform the LCFL system into a **left-linear** system over a **larger** semiring (the tensor product semiring). Solve with Tarjan. Transform back.

The tensor product trick:

```
a ⊗ Y ⊗ b  →  Z ⊗_T (a^t ⊗ b)
```

Now Z is only on the LEFT! The sandwich has been "uncurried" into the tensor structure.

## Why Newton + Tensor Works

1. **Newton** reduces the nonlinear problem to a sequence of linear problems
2. **Tensor products** reduce each LCFL linear problem to a regular problem
3. **Tarjan** solves each regular problem efficiently
4. Repeat for ~3-4 rounds until convergence

Each round is polynomial work. Number of rounds is constant. Total: polynomial time for what would otherwise be exponential.

## Summary

| Problem | Solution |
|---------|----------|
| Loops within a procedure | Kleene star in path expressions |
| Non-regular path structure | Tarjan's algorithm |
| Interprocedural sandwiches | Newton's method (linearize) |
| Linearized system still LCFL | Tensor products (regularize) |

The full algorithm (NPA-TP) combines all of these. It's beautiful machinery.

## Problem Set 6

### Problem 6.1: Sandwich Equations

Consider:
```
X = a ⊗ X ⊗ b ⊕ c
```

With `a`, `b`, `c` being 2×2 Boolean matrices:
- `a = [[1,0],[0,0]]` (projects to first component)
- `b = [[0,0],[0,1]]` (projects to second component)  
- `c = [[0,1],[0,0]]` (maps second to first)

What is `a ⊗ c ⊗ b`? Is it zero?

Now compute a few iterations of Kleene starting from `X = 0`:
- X^(0) = 0
- X^(1) = f(X^(0)) = c
- X^(2) = f(X^(1)) = a ⊗ c ⊗ b ⊕ c = ?

Does the sequence stabilize?

### Problem 6.2: Derivatives

Compute the formal derivative of `X² = X ⊗ X` with respect to X.

Remember: product rule gives `D[X ⊗ X] = D[X] ⊗ X ⊕ X ⊗ D[X] = Y ⊗ X ⊕ X ⊗ Y`

If X is a 2×2 matrix and the derivative is with respect to X, what's the "type" of Y ⊗ X ⊕ X ⊗ Y? (Hint: it's still a 2×2 matrix, but it's linear in Y)

### Problem 6.3: Star Derivative

For the expression `X*`, the derivative is `X* ⊗ Y ⊗ X*`.

Verify this makes sense by expanding `X* = 1 ⊕ X ⊕ X² ⊕ X³ ⊕ ...` and differentiating term by term:
- D[1] = 0
- D[X] = Y
- D[X²] = Y ⊗ X ⊕ X ⊗ Y
- D[X³] = ?

Sum these up. Do you see the pattern `X* ⊗ Y ⊗ X*` emerging?

### Problem 6.4: Why LCFL?

A **context-free grammar** generates strings using rules like `S → aSb | c`.

An **LCFL** (Linear Context-Free Language) is generated by a grammar where each rule has at most one non-terminal:
- `S → aSb | c` (linear — one S per rule)
- `S → aSbSc` (nonlinear — two S's)

The equation `X = a ⊗ X ⊗ b ⊕ c` corresponds to `S → aSb | c`.

What language does this generate? Write out a few strings.

### Problem 6.5: Left-Linear vs LCFL

Compare these two equation systems:

Left-linear: `X = a ⊗ X ⊕ c`
LCFL: `X = a ⊗ X ⊗ b ⊕ c`

For the left-linear system, we can solve directly: `X = a* ⊗ c`.

Verify this is correct by substituting back:
```
a* ⊗ c = a ⊗ (a* ⊗ c) ⊕ c  ?
```

(Hint: use `a* = 1 ⊕ a ⊗ a*`)

Why doesn't the same trick work for LCFL?

### Problem 6.6: Mutual Recursion

Consider:
```
X = Y ⊗ a
Y = X ⊗ b ⊕ c
```

Substitute Y into X:
```
X = (X ⊗ b ⊕ c) ⊗ a = X ⊗ b ⊗ a ⊕ c ⊗ a
```

This is RIGHT-linear in X. Solve it.

Now consider:
```
X = a ⊗ Y
Y = b ⊗ X ⊕ c
```

Substitute Y into X:
```
X = a ⊗ (b ⊗ X ⊕ c) = a ⊗ b ⊗ X ⊕ a ⊗ c
```

This is LEFT-linear in X. Solve it.

What about:
```
X = a ⊗ Y ⊗ b
Y = c ⊗ X ⊗ d ⊕ e
```

Substitute. What form does the equation for X take?

### Problem 6.7: Iteration Count

For left-linear `X = a ⊗ X ⊕ c` with Boolean matrices, how many iterations until convergence?

For sandwich `X = a ⊗ X ⊗ b ⊕ c`, the iteration can take longer. Can you construct an example where it takes more than 2 iterations to converge?

### Problem 6.8: Newton Intuition

In ordinary calculus, Newton's method for f(x) = x² - 2 (finding √2) converges quadratically: each iteration doubles the correct digits.

For semiring equations, convergence is also fast but for a different reason. The "derivative" captures the structure of one layer of recursion. After ~log(depth) rounds, you've captured all the recursion.

Think about `X = a ⊗ X ⊗ b ⊕ c` as building a "tree" of nested contexts. How does the Newton iteration "peel off" layers of this tree?

---

Next: [Lecture 7: Going Interprocedural with Tensor Products](./07-tensor-products.md)
