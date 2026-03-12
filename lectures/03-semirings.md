# Lecture 3: Semirings Are Everywhere

Last lecture we built up fixed points and lattices. Now let's take a step back and notice something remarkable: a bunch of seemingly different problems all have the same algebraic structure.

## Three Problems

**Problem 1: Shortest Paths**

Given a weighted graph, find the shortest path from source to each node.

At each node, we track the minimum distance from source. When paths merge at a node, we take the **minimum**. When we traverse an edge, we **add** the weight.

**Problem 2: Reachability**

Given a graph, determine which nodes are reachable from a source.

At each node, we track "is this reachable?" (a boolean). When paths merge, we take **or**. When we traverse an edge, we take **and** with the edge's existence.

**Problem 3: Reaching Definitions**

At each node, track which definitions might reach here.

When paths merge, we take **union**. When we traverse an edge, we apply a transfer function (gen/kill).

## The Pattern

Notice anything similar?

| Problem | "Addition" | "Multiplication" | "Zero" | "One" |
|---------|------------|-----------------|--------|-------|
| Shortest paths | min | + | ∞ | 0 |
| Reachability | or | and | false | true |
| Reaching defs | ∪ | (transfer function) | ∅ | (identity function) |

In each case:
- There's a way to **combine** information from multiple paths (the "addition")
- There's a way to **extend** information along a path (the "multiplication")
- There's an identity for combining (the "zero" — contributes nothing)
- There's an identity for extending (the "one" — doesn't change anything)

This is a **semiring**.

## Definition: Semiring

A semiring (S, ⊕, ⊗, 0, 1) is a set S with two operations:

**Combine** (⊕):
- Associative: (a ⊕ b) ⊕ c = a ⊕ (b ⊕ c)
- Commutative: a ⊕ b = b ⊕ a
- Has identity 0: a ⊕ 0 = a

**Extend** (⊗):
- Associative: (a ⊗ b) ⊗ c = a ⊗ (b ⊗ c)
- Has identity 1: a ⊗ 1 = 1 ⊗ a = a
- Zero annihilates: a ⊗ 0 = 0 ⊗ a = 0

**Connection**:
- Extend distributes over combine: a ⊗ (b ⊕ c) = (a ⊗ b) ⊕ (a ⊗ c)
- And on the other side: (b ⊕ c) ⊗ a = (b ⊗ a) ⊕ (c ⊗ a)

If this looks like a ring without subtraction, that's exactly what it is. The "semi" means we're missing additive inverses.

## Instantiating the Semirings

**Shortest path semiring** (min, +):
- S = ℕ ∪ {∞} (non-negative integers plus infinity)
- ⊕ = min
- ⊗ = +
- 0 = ∞ (minimum with ∞ gives the other thing)
- 1 = 0 (adding 0 doesn't change the distance)

**Boolean semiring** (or, and):
- S = {false, true}
- ⊕ = or
- ⊗ = and
- 0 = false
- 1 = true

**Powerset semiring** for sets of paths:
- S = 2^Paths (sets of paths)
- ⊕ = union
- ⊗ = path concatenation (pairwise)
- 0 = ∅ (no paths)
- 1 = {ε} (set containing just the empty path)

## Why This Abstraction Helps

Once you see the pattern, you can write algorithms that work for ANY semiring. The dataflow algorithm becomes:

```
for each node n:
    X[n] := ⊕ { edge_weight(p, n) ⊗ X[p] | p predecessor of n }
```

Same algorithm, different semirings = different analyses.

## Kleene Star

There's one more operation we need: the **Kleene star**.

a* = 1 ⊕ a ⊕ (a ⊗ a) ⊕ (a ⊗ a ⊗ a) ⊕ ...

It's the infinite sum of all powers of a. This represents "going around a loop zero or more times."

For different semirings:
- Boolean: true* = true, false* = true
- (min, +): 0* = 0, anything else* = 0 (zero iterations is best)
- Sets of paths: a* = {ε, a, aa, aaa, ...}

The star is what lets us handle loops without infinite iteration.

## Idempotent Semirings

A semiring is **idempotent** if a ⊕ a = a for all a.

Boolean is idempotent: true or true = true.
(min, +) is idempotent: min(5, 5) = 5.
But ({1,2}, {1,2,3}) under union is idempotent too.

Idempotence matters because it means information stabilizes. Once you've seen a definition, seeing it again doesn't change anything. This is why dataflow iteration terminates.

## Example: Graph Reachability

Let's work through reachability concretely.

Graph:
```
1 → 2 → 3
↓   ↑
4 → 5
```

Question: Which nodes can reach node 3?

Using the Boolean semiring (or, and):
- Initialize: reach[3] = true (node 3 reaches itself)
- Propagate backwards (we're computing "can reach", not "is reachable from")

After propagation:
- Node 2 can reach 3: true (direct edge)
- Node 1 can reach 3: true (via 2)
- Node 4 can reach 3: true (via 5 → 2 → 3)
- Node 5 can reach 3: true (via 2)

The semiring operations:
- ⊕ = or: if ANY predecessor can reach the target, we can too
- ⊗ = and: we can reach target if the edge exists AND our successor can reach it

## Example: Shortest Paths

Same graph with weights:
```
1 →(2)→ 2 →(1)→ 3
↓(5)    ↑(3)
4 →(1)→ 5
```

Using (min, +):
- Initialize: dist[1] = 0 (source)
- Propagate forward

After iteration:
- dist[2] = min(0+2, 0+5+1+3) = min(2, 9) = 2
- dist[3] = 2+1 = 3
- dist[4] = 0+5 = 5
- dist[5] = 5+1 = 6

The semiring operations:
- ⊕ = min: take the best path
- ⊗ = +: accumulate distances

## Matrix Semirings

Here's where it gets interesting for dataflow.

An N×N matrix over a semiring S is itself a semiring:
- ⊕ = pointwise combine (OR for boolean matrices)
- ⊗ = matrix multiplication (using S's operations)
- 0 = zero matrix
- 1 = identity matrix

For Boolean matrices, ⊗ is relational composition: (A ⊗ B)[i,j] = ∃k. A[i,k] ∧ B[k,j]

This is huge! It means we can represent **relations** between program states as matrices, and compose them using matrix multiplication.

## Preview: Predicate Abstraction

In predicate abstraction, we track which predicates hold at each program point.

With k predicates, there are 2^k possible predicate assignments. A program statement relates pre-states to post-states — that's a 2^k × 2^k Boolean matrix!

- Entry to node: which predicates might hold?
- Edge effect: how does the statement change predicates?
- Composition: matrix multiplication!

The matrix semiring is the foundation of interprocedural analysis. We'll return to this in later lectures.

## Why Not Just Use Sets?

You might ask: why this semiring abstraction instead of just "sets of facts"?

Because extend (⊗) isn't always union or intersection. For relational analyses, ⊗ is composition. For weighted analyses, ⊗ accumulates costs. The semiring abstraction captures the essential structure without committing to a specific domain.

This lets us write algorithms once and instantiate them for different analyses.

## Connection to Last Lecture

Remember lattices from Lecture 2? Semirings and lattices are related:

- ⊕ induces a partial order: a ≤ b iff a ⊕ b = b
- For idempotent semirings, this makes S a join-semilattice with ⊕ as join
- The fixed-point machinery from Lecture 2 applies

But semirings give us MORE than lattices: the extend operation (⊗) that models path composition. This is what enables the path-expression approach we're building toward.

## In Code

Here's how we define semirings in Rust:

```rust
pub trait Semiring: Clone + Eq {
    fn zero() -> Self;
    fn one() -> Self;
    fn combine(&self, other: &Self) -> Self;  // ⊕
    fn extend(&self, other: &Self) -> Self;   // ⊗
    fn star(&self) -> Self;                   // Kleene star
}
```

And the Boolean semiring:

```rust
impl Semiring for bool {
    fn zero() -> Self { false }
    fn one() -> Self { true }
    fn combine(&self, other: &Self) -> Self { *self || *other }
    fn extend(&self, other: &Self) -> Self { *self && *other }
    fn star(&self) -> Self { true }  // a* = 1 ⊕ a ⊕ ... = true for any boolean
}
```

## Problem Set 3

### Problem 3.1: Verify the Axioms

Check that the (min, +) semiring actually satisfies all the semiring axioms:

a) min is associative, commutative, with identity ∞
b) + is associative with identity 0
c) min(a, ∞) = a and a + ∞ = ∞ (zero annihilates)
d) + distributes over min: a + min(b, c) = min(a+b, a+c)

### Problem 3.2: A Broken Semiring

Someone proposes the "maximum-product" semiring for graphs:
- ⊕ = max
- ⊗ = ×
- 0 = 0
- 1 = 1

Is this actually a semiring? Check the axioms. (Hint: what's the identity for max?)

### Problem 3.3: String Semiring

Define a semiring over strings where:
- ⊕ = shortest string (or lexicographic for ties)
- ⊗ = concatenation

What's the identity for concatenation? What string should be 0? Does this work?

### Problem 3.4: Matrix Verification

For 2×2 Boolean matrices, verify that:
- The zero matrix is the identity for ⊕ (pointwise OR)
- The identity matrix is the identity for ⊗ (matrix multiply)
- Zero annihilates: M ⊗ 0 = 0 ⊗ M = 0

### Problem 3.5: Kleene Star for Matrices

For a 2×2 Boolean matrix A, the star A* should represent "zero or more applications of A."

Given A = [[0,1],[1,0]] (swap), compute A⁰, A¹, A², A³. What's the pattern?

Now compute A* = A⁰ ⊕ A¹ ⊕ A² ⊕ ...

### Problem 3.6: Weighted Paths

Consider this weighted graph:
```
    2
A ───→ B
│       │
3       1
↓       ↓
C ───→ D
    4
```

Using the (min, +) semiring with A as source:
a) What's the shortest path to each node?
b) Write this as matrix operations (4×4 matrix, one entry per edge weight, run powers until stable)

### Problem 3.7: Non-Commutative Extend

The Boolean semiring has commutative ⊗ (and is commutative). But matrix multiplication is NOT commutative.

Give two 2×2 Boolean matrices A and B where A ⊗ B ≠ B ⊗ A.

Why does this non-commutativity matter for dataflow? (Hint: think about the order of edges on a path)

### Problem 3.8: Relating Lattices and Semirings

For an idempotent semiring S, the relation a ≤ b ⟺ a ⊕ b = b is a partial order.

a) Verify this for the Boolean semiring
b) Verify this for the (min, +) semiring
c) Show that ⊕ is the join (least upper bound) in this order

---

Next: [Lecture 4: Path Expressions](./04-path-expressions.md)
