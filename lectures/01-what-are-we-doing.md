# Lecture 1: What Are We Even Doing?

Before we write any code, we need to understand the problem we're trying to solve. This lecture sets up the motivation for everything that follows.

## The Setup

You have a program. You want to know things about it *without running it*. Things like:

- "Which variable definitions might reach this point?"
- "Is this variable always initialized before use?"
- "Could this pointer be null here?"
- "What values might this expression have?"

This is **static analysis** — analyzing programs at compile time rather than runtime.

## A Concrete Example

Consider this tiny program:

```
1: x = 5
2: y = x + 1
3: if (condition) {
4:   x = 10
5: }
6: z = x + y
```

**Question**: At line 6, which definitions of `x` might be "live"? In other words, which assignments to `x` could have been the most recent one?

Think about it. There are two paths through this program:
- Path A: 1 → 2 → 3 → 6 (skip the if-body)
- Path B: 1 → 2 → 3 → 4 → 5 → 6 (take the if-body)

On path A, the definition at line 1 reaches line 6.
On path B, the definition at line 4 reaches line 6.

So the answer is: **both definitions might reach line 6**. We don't know which path will be taken at runtime, so we have to be conservative and say "either one."

This is called **reaching definitions** analysis.

## The Control Flow Graph

We represent programs as **control flow graphs** (CFGs). Each node is a program point, each edge is a possible transition. Our example becomes:

```
    [1: x=5]
        |
        v
    [2: y=x+1]
        |
        v
    [3: if]
      /    \
     v      v
[4: x=10]  [skip]
     |      |
     v      v
    [5: endif]
        |
        v
    [6: z=x+y]
```

(I'm simplifying — the exact structure depends on how you model conditionals. The point is: there's a graph, and we care about paths through it.)

## What We're Computing

At each node, we want to compute some **information** based on all possible paths that could reach it. For reaching definitions:

- At node 1: the definition `x=5` is created here
- At node 2: the definition `x=5` reaches here (from node 1)
- At node 3: the definition `x=5` still reaches here
- At node 4: the definition `x=10` is created here, kills `x=5`
- At node 5: either `x=5` or `x=10` might reach here (depending on path)
- At node 6: either `x=5` or `x=10` might reach here

The key insight: **we're combining information from all paths**.

## The Fixed-Point View

Here's the standard way to think about this. For each node n, define:

```
IN[n] = information flowing into n (from predecessors)
OUT[n] = information flowing out of n (after n's effects)
```

For reaching definitions:
- OUT[n] = (IN[n] - killed definitions) ∪ (definitions created at n)
- IN[n] = ∪ { OUT[p] | p is a predecessor of n }

The second equation says: the information flowing into n is the **union** of information from all predecessors.

Now here's the problem: these equations are **circular**. To compute IN[n], you need OUT of predecessors. To compute OUT of predecessors, you need their IN. And if there are loops, you're computing things in terms of themselves.

## Iteration to the Rescue

The classic solution: **iterate until nothing changes**.

```
initialize all IN[n] and OUT[n] to empty
repeat:
    for each node n:
        IN[n] = ∪ { OUT[p] | p predecessor of n }
        OUT[n] = transfer(IN[n])
until nothing changed
```

This is called **Kleene iteration** or **chaotic iteration**. It works! For many analyses, it converges to the right answer.

But there's a problem lurking.

## The Problem With Nested Loops

Consider a program with nested loops:

```
while (...) {
    while (...) {
        while (...) {
            x = ...
        }
    }
}
```

Each loop level can multiply the number of iterations needed. If you have k levels of nesting, you might need O(n^k) iterations in the worst case.

For a simple reachability analysis on a path graph with a single back-edge, you need O(n) iterations — one for each node. Add another loop inside, and you might need O(n²). And so on.

This is not a theoretical concern. Real programs have deeply nested loops. Analyzing them with naive iteration can be *slow*.

## A Different Approach

What if instead of iterating, we could just... **write down the answer**?

Think about it differently. The information at node n depends on all paths from entry to n. Each path transforms the information in some way. We want the combination of all those transformations.

If we could express "all paths from entry to n" as a **formula**, and then **evaluate** that formula, we'd have our answer directly. No iteration needed.

This is the path-expression approach. We'll build up to it over the next few lectures.

## What's Coming

- **Lecture 2**: We'll make "combine information from paths" precise using lattices and fixed points.

- **Lecture 3**: We'll notice that different analyses (shortest paths, reachability, dataflow) all have the same algebraic structure. This is the semiring abstraction.

- **Lecture 4**: We'll learn to write down "all paths" as regular expressions.

- **Lecture 5**: We'll use Tarjan's algorithm to compute path expressions efficiently.

- **Lecture 6**: We'll see why interprocedural analysis breaks path expressions — and introduce Newton's method as the fix.

- **Lecture 7**: We'll handle function calls with tensor products.

- **Lecture 8**: We'll implement the whole thing in Rust.

## Problem Set 1

These problems don't require any code yet. They're about building intuition.

### Problem 1.1: More Reaching Definitions

Consider this program:

```
1: x = 1
2: y = 2
3: while (condition) {
4:     x = y
5:     y = x
6: }
7: z = x + y
```

Which definitions of `x` might reach line 7? Which definitions of `y`? Trace through a few iterations of the loop to convince yourself.

### Problem 1.2: A Different Analysis

**Available expressions** asks: which expressions have *definitely* been computed on *all* paths to this point (and not invalidated)?

```
1: a = b + c
2: if (condition) {
3:     d = b + c
4: } else {
5:     e = f + g
6: }
7: x = b + c
```

At line 7, is the expression `b + c` available? What about `f + g`?

How does this analysis differ from reaching definitions? (Hint: "might reach" vs "definitely available")

### Problem 1.3: The Cost of Iteration

Consider a simple path graph: 1 → 2 → 3 → ... → n, with a back-edge from n to 1.

If we're computing reachability (can node 1 reach node k?), how many iterations of Kleene iteration do we need before the information "node 1 is reachable" propagates to node n?

What if we add another back-edge from n/2 to 1? Does this change the number of iterations?

### Problem 1.4: Thinking About Structure

In the loop example:
```
while (...) {
    while (...) {
        ...
    }
}
```

The outer loop's body contains the entire inner loop. When we iterate, information from the inner loop has to "escape" to the outer loop, which then has to iterate again.

Can you think of a program structure where the iteration count grows *exponentially* with the nesting depth? (This is tricky — don't worry if you can't get it, we'll cover it later.)

### Problem 1.5: Why Undecidability Doesn't Stop Us

You might know that many program properties are undecidable — we can't always give the exact right answer. Yet dataflow analysis "works." How?

Think about the difference between:
- "Definition D definitely reaches point P"
- "Definition D might reach point P"
- "Definition D definitely does NOT reach point P"

Which of these can we compute exactly? Which do we approximate, and in which direction?

---

Next: [Lecture 2: Fixed Points and Iteration](./02-fixed-points.md)
