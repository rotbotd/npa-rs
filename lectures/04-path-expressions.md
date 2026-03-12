# Lecture 4: Path Expressions

We've seen that dataflow analysis is about computing information over all paths from entry to each node. This lecture introduces path expressions — regular expressions that describe all those paths.

## The Key Idea

Instead of *iterating* to find the fixed point, what if we could *write down* a formula for "all paths from entry to node N"?

Consider this CFG:
```
1 → 2 → 3
```

The paths from 1 to 3: just `1 → 2 → 3`. We can write this as `a ⊗ b` where `a` labels the first edge and `b` labels the second.

Now add a loop:
```
1 → 2 → 3 → 2 (back to 2)
    ↓
    4
```

The paths from 1 to 4:
- `1 → 2 → 4` (skip the loop)
- `1 → 2 → 3 → 2 → 4` (one iteration)
- `1 → 2 → 3 → 2 → 3 → 2 → 4` (two iterations)
- ...

Labeling edges: `a` (1→2), `b` (2→3), `c` (3→2), `d` (2→4)

All paths to 4: `a ⊗ d ⊕ a ⊗ b ⊗ c ⊗ d ⊕ a ⊗ b ⊗ c ⊗ b ⊗ c ⊗ d ⊕ ...`

This simplifies to: `a ⊗ (b ⊗ c)* ⊗ d`

That's a **path expression** — a regular expression over edge labels.

## Regular Expressions Over Semirings

Path expressions are regular expressions, but instead of string matching, we evaluate them in a semiring.

**Syntax**:
```
e ::= 0           (zero — no paths)
    | 1           (one — empty path)
    | a           (constant — single edge)
    | e ⊕ e       (combine — alternative paths)
    | e ⊗ e       (extend — sequential composition)
    | e*          (star — zero or more repetitions)
```

**Evaluation** in a semiring S:
- `eval(0) = S.zero`
- `eval(1) = S.one`
- `eval(a) = a` (where a ∈ S)
- `eval(e₁ ⊕ e₂) = eval(e₁) ⊕ eval(e₂)`
- `eval(e₁ ⊗ e₂) = eval(e₁) ⊗ eval(e₂)`
- `eval(e*) = eval(e)*`

The beautiful thing: same expression, different semirings, different meanings.

## Examples

**Boolean semiring** (reachability):
- Edge labels: all `true`
- `a ⊗ (b ⊗ c)* ⊗ d` evaluates to `true ∧ true* ∧ true = true`
- Meaning: node 4 is reachable from node 1

**(min, +) semiring** (shortest paths):
- Edge labels: weights `a=2, b=1, c=1, d=3`
- `a ⊗ (b ⊗ c)* ⊗ d = 2 + (1+1)* + 3`
- Star of 2 is 0 (zero iterations is best)
- Result: `2 + 0 + 3 = 5`

**Matrix semiring** (predicate abstraction):
- Edge labels: transition relations as N×N Boolean matrices
- `a ⊗ (b ⊗ c)* ⊗ d` computes the relation from entry predicates to exit predicates
- Star = transitive closure of the loop body

## Why This Helps

**Kleene iteration**: Run the loop 0 times, then 1 time, then 2 times... hope it converges.

**Path expressions**: Write `(bc)*` and compute the star *directly*. For Boolean matrices, this is transitive closure — O(n³) once, not O(n³ × iterations).

The win is especially big when:
1. Star is cheaper than repeated iteration (matrices, regular languages)
2. There are nested loops (composition beats nesting)

## The Catch

Computing path expressions naively is expensive. For a general graph with n nodes and m edges, Kleene's algorithm (the classic RE construction) is O(n³).

And worse: the expressions can grow exponentially in size if you're not careful about sharing.

## Preview: Tarjan's Algorithm

Next lecture covers Tarjan's path-expression algorithm, which exploits **dominators** to compute path expressions efficiently.

The key insight: if every path to node N goes through node D, we can decompose:
```
path(N) = path(D) ⊗ local_path(D → N)
```

This recursive decomposition yields nearly-linear time for reducible CFGs.

## From Paths to Dataflow

Once we have path expressions, dataflow analysis becomes:

1. Compute path expression for each node
2. Evaluate in the appropriate semiring

For reaching definitions:
- Edge labels: transfer functions (what definitions are gen'd/kill'd)
- Compose via function composition (⊗)
- Combine via union (⊕)
- Star: fixed point of loop body

The path expression directly encodes the dataflow equations!

## Variables in Path Expressions

For interprocedural analysis, path expressions contain *variables* representing procedure summaries:

```
X₁ = a ⊗ X₂ ⊗ b ⊕ c
X₂ = d ⊗ X₁ ⊗ e ⊕ f
```

Here `X₁` and `X₂` are the summaries for two mutually recursive procedures.

These are **polynomial equations** over the semiring. Solving them is harder than just evaluating a star — that's where Newton's method comes in (Lecture 6).

## Problem Set 4

### Problem 4.1: Path Expression by Hand

For this CFG:
```
    1
   / \
  a   b
  ↓   ↓
  2   3
   \ /
    c d
    ↓ ↓
    4
```
(edges a: 1→2, b: 1→3, c: 2→4, d: 3→4)

Write the path expression for paths from 1 to 4.

### Problem 4.2: Loop Path Expression

For this CFG with nested loops:
```
1 → 2 → 3 → 4 → 3 (inner loop)
        ↓
        5 → 2 (outer loop)
        ↓
        6
```
Edges: a(1→2), b(2→3), c(3→4), d(4→3), e(3→5), f(5→2), g(5→6)

Write path expressions for nodes 3, 5, and 6.

### Problem 4.3: Evaluation in Boolean Semiring

Take your answer from 4.1 and evaluate it in the Boolean semiring (where all edges are `true`).

What's the result? What does it mean?

### Problem 4.4: Evaluation in (min, +)

Same CFG as 4.1, but now with weights: a=1, b=2, c=3, d=1.

Evaluate the path expression. What's the shortest path to node 4?

### Problem 4.5: Expression Size

Consider a "ladder" graph:
```
1 → 2 → 3 → ... → n
↓   ↓   ↓         ↓
1'→ 2'→ 3'→ ... → n'
```

With horizontal edges aᵢ, bᵢ and vertical edges cᵢ.

How many distinct paths are there from 1 to n'? (Think about this combinatorially.)

If we write out the path expression without any sharing, how big would it be?

This illustrates why naive path expression construction blows up.

### Problem 4.6: Star Semantics

For the (min, +) semiring:
- What is 0* (star of the additive identity)?
- What is 1* (star of the multiplicative identity)?
- What is 5* (star of a positive number)?

Justify your answers using the definition a* = 1 ⊕ a ⊕ a² ⊕ a³ ⊕ ...

### Problem 4.7: Matrix Path Expression

Consider a 2×2 Boolean matrix semiring.

Let A = [[1,0],[1,1]] (identity plus edge 0→1).

Compute A*, A², A³. What pattern emerges?

What does A* represent in terms of reachability?

### Problem 4.8: The Interprocedural Challenge

Consider these equations:
```
X = a ⊗ Y ⊗ b
Y = c ⊗ X ⊗ d ⊕ e
```

Try to "unfold" X by substituting the equation for Y:
```
X = a ⊗ (c ⊗ X ⊗ d ⊕ e) ⊗ b
```

Can you simplify this to a form `X = something ⊗ X ⊗ something_else ⊕ constant`?

Notice how X appears "sandwiched" — this is the non-commutativity problem we'll tackle in Lecture 7.

---

Next: [Lecture 5: Tarjan's Algorithm](./05-tarjan.md)
