# Lecture 2: Fixed Points and Iteration

Last time we saw that dataflow analysis involves circular equations: the information at a node depends on its predecessors, which might depend on it. This lecture makes that precise.

## The Equation System

Let's write down what we're computing more formally. For each node n in the CFG:

```
X[n] = F_n(X[pred₁], X[pred₂], ..., X[predₖ])
```

Where:
- `X[n]` is the dataflow information at node n
- `pred₁, ..., predₖ` are the predecessors of n
- `F_n` combines the information from predecessors and applies n's transfer function

For reaching definitions, `F_n` looks like:

```
F_n(in₁, in₂, ..., inₖ) = GEN[n] ∪ ((in₁ ∪ in₂ ∪ ... ∪ inₖ) - KILL[n])
```

Where:
- `GEN[n]` = definitions created at n
- `KILL[n]` = definitions killed at n

## What Does "Solution" Mean?

We want values `X[n]` for every node such that all equations are satisfied simultaneously. Such a solution is called a **fixed point** — if you plug the values back into the equations, you get the same values out.

But here's the rub: there might be **multiple fixed points**!

Consider a trivial example:

```
X = X
```

Any value of X is a fixed point. That's not helpful.

Or consider:

```
X = X ∪ {a}
```

The only fixed point is a set containing `a`. But how big? `{a}`, `{a, b}`, and `{a, b, c, ...}` are all fixed points.

## We Want the Smallest One

For dataflow analysis, we want the **least fixed point** — the smallest solution that satisfies all equations.

Why? Because we're approximating. If we say "definition d might reach point p," we want that to be true only if there's actually a path that could make it happen. Extra information is spurious.

The least fixed point gives us the tightest correct approximation.

## Enter Lattices

To make "smallest" precise, we need an ordering on our domain. This is where lattices come in.

A **partially ordered set** (poset) is a set with a relation ≤ that is:
- Reflexive: x ≤ x
- Antisymmetric: if x ≤ y and y ≤ x, then x = y
- Transitive: if x ≤ y and y ≤ z, then x ≤ z

A **lattice** is a poset where every two elements have:
- A **least upper bound** (join, ⊔): smallest element ≥ both
- A **greatest lower bound** (meet, ⊓): largest element ≤ both

A **complete lattice** has joins and meets for all subsets, including infinite ones. It has:
- A **bottom** element ⊥: smallest element
- A **top** element ⊤: largest element

## Example: Sets Under Inclusion

The powerset of some universe U, ordered by ⊆, is a complete lattice:
- ⊔ = union
- ⊓ = intersection
- ⊥ = empty set
- ⊤ = U

For reaching definitions with definitions {d₁, d₂, d₃}, the lattice is:

```
         {d₁, d₂, d₃}          ← top (all definitions)
        /     |      \
   {d₁,d₂}  {d₁,d₃}  {d₂,d₃}
      |    \  /  \   /   |
     {d₁}   {d₂}   {d₃}
        \    |    /
          empty              ← bottom (no definitions)
```

## Monotonicity

For the least fixed point to exist and be computable, our functions need to be **monotonic**:

If x ≤ y, then f(x) ≤ f(y)

In other words: more input gives at least as much output. Functions that satisfy this "play nice" with the ordering.

Dataflow transfer functions are almost always monotonic. If more definitions flow into a node, at least as many flow out.

## The Fixed Point Theorem

**Kleene's Fixed Point Theorem**: If L is a complete lattice and f: L → L is monotonic, then:

1. f has a least fixed point
2. It equals ⊔{fⁿ(⊥) | n ≥ 0} = ⊥ ⊔ f(⊥) ⊔ f(f(⊥)) ⊔ ...

In other words: start at bottom, keep applying f, take the limit. You'll reach the least fixed point.

This is the theoretical justification for Kleene iteration!

## Iteration in Practice

```python
def kleene_iteration(equations, bottom):
    X = [bottom for _ in equations]  # start at ⊥
    changed = True
    while changed:
        changed = False
        for i, eq in enumerate(equations):
            new_val = eq.evaluate(X)
            if new_val != X[i]:
                X[i] = new_val
                changed = True
    return X
```

Each iteration computes f(f(...f(⊥)...)). We stop when nothing changes — that's the fixed point.

## When Does It Terminate?

If the lattice has **finite height** (longest ascending chain is finite), iteration always terminates. Each step either stays the same or goes up in the lattice. Can't go up forever.

For reaching definitions on a finite program, the lattice has height = number of definitions. We're guaranteed to terminate in at most that many iterations per node.

But what if the lattice has infinite height? We'll need **widening** — a topic for later.

## Why Multiple Iterations?

Consider:

```
1: x = 1
2: y = x
3: if (...) goto 2
4: z = y
```

CFG:
```
1 → 2 → 3 → 4
    ↑   |
    +---+
```

On first iteration:
- At node 1: {x=1}
- At node 2: {x=1} (from node 1)... but wait, what about from node 3?

We don't know what reaches node 3 yet! We haven't computed it. So on first pass, we might underestimate.

On second iteration:
- At node 2: {x=1} from node 1, plus whatever comes from node 3
- At node 3: information from node 2 (now has more)
- At node 4: information from node 3 (now has more)

Information "propagates" through the graph. Each iteration pushes it one edge further.

## The Cost of Propagation

Here's the key insight that motivates everything to come:

In a straight-line program (no loops), information propagates in one pass. You just process nodes in topological order.

But loops require multiple passes. Information has to go around the loop, then around again, then again...

**Nested loops are worse.** The inner loop might need many iterations. Each iteration of the outer loop restarts the inner loop. Multiply the costs.

For k levels of loop nesting with n nodes each, naive iteration can take O(nᵏ) time.

This is why we'll eventually want something smarter than iteration.

## Connecting to Last Lecture

Remember Problem 1.4 asked about exponential blowup? Here's the answer:

```
while (...) {           // loop 1
    while (...) {       // loop 2
        while (...) {   // loop 3
            x = ...
        }
    }
}
```

Each level of nesting can multiply the iteration count. Three loops deep with 10 nodes each: potentially 10³ = 1000 iterations.

The theory says it terminates. It doesn't say it terminates quickly.

## Preview: Path Expressions

Next lecture we'll introduce semirings — the algebraic abstraction that unifies different dataflow analyses.

Then we'll learn to write down "all paths" as a formula. Instead of iterating, we'll compute the answer directly. This is the path-expression approach.

For now, the key takeaway: **we're computing least fixed points over lattices**, and iteration works but can be slow.

## Problem Set 2

### Problem 2.1: Drawing Lattices

Draw the lattice for the powerset of {a, b, c} ordered by ⊆. Label all 8 elements and draw edges for the covering relation (x covers y if x > y and there's no z with x > z > y).

### Problem 2.2: Monotonicity Check

Which of these functions on sets are monotonic (where ≤ means ⊆)?

a) f(X) = X ∪ {a}
b) f(X) = X - {a}  
c) f(X) = {a} if |X| > 3, else ∅
d) f(X) = X ∩ Y for some fixed Y

For each non-monotonic one, give a counterexample.

### Problem 2.3: Computing Fixed Points

Consider the equations over sets:

```
X = {a} ∪ Y
Y = X ∩ {a, b}
```

Starting from X = Y = ∅, trace Kleene iteration until convergence. How many iterations does it take?

What is the least fixed point?

### Problem 2.4: Iteration Order

In the iteration algorithm, we process equations in some order. Does the order affect:

a) Whether we find the correct fixed point?
b) How many iterations we need?

Hint: Consider a chain A → B → C → D. What if we process in order D, C, B, A versus A, B, C, D?

### Problem 2.5: Infinite Lattices

Consider the lattice of natural numbers with ≤ (so 0 ≤ 1 ≤ 2 ≤ ...).

a) Does this lattice have a top element?
b) Is it a complete lattice? (Hint: consider the join of all natural numbers)
c) Consider f(n) = n + 1. Is it monotonic? Does it have a fixed point?

This is why we usually need **finite height** or special techniques for infinite domains.

### Problem 2.6: Available Expressions Lattice

Available expressions analysis asks: which expressions are definitely computed on ALL paths?

a) Should the lattice order be "more expressions = higher" or "more expressions = lower"?
b) What should the initial value at the entry node be?
c) How do you combine information from multiple predecessors — union or intersection?

Hint: "definitely on all paths" means we need something to be true on EVERY predecessor.

### Problem 2.7: The "May" vs "Must" Distinction

Reaching definitions is a "may" analysis: definition d MAY reach point p.
Available expressions is a "must" analysis: expression e MUST be available at point p.

How does this distinction affect:
- The direction of the lattice ordering
- The choice of ⊔ vs ⊓ for combining paths
- Whether we want the least or greatest fixed point

---

Next: [Lecture 3: Semirings Are Everywhere](./03-semirings.md)
