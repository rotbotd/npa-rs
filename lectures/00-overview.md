# From Zero to Newtonian Program Analysis

A series of lectures and problem sets that take you from "I know what a control flow graph is" to implementing a full interprocedural dataflow analysis engine.

## Prerequisites

- You know what a control flow graph (CFG) is
- You can read and write basic Rust (or are willing to learn as you go)
- You're comfortable with sets, functions, and basic discrete math

That's it. We'll build everything else from scratch.

## The Journey

1. [**What Are We Even Doing?**](./01-what-are-we-doing.md) — The dataflow analysis problem, why we care, what "reaching definitions" means
2. [**Fixed Points and Iteration**](./02-fixed-points.md) — Discovering lattices by accident, iterating until things stop changing
3. [**Semirings Are Everywhere**](./03-semirings.md) — Abstracting the pattern: shortest paths, reachability, dataflow are all the same thing
4. [**Path Expressions**](./04-path-expressions.md) — From "run the algorithm" to "write down the answer" 
5. [**Tarjan's Algorithm**](./05-tarjan.md) — Computing path expressions efficiently with dominators
6. [**The Interprocedural Problem**](./06-interprocedural.md) — Why sandwiches break everything, and Newton's method as the fix
7. [**Going Interprocedural with Tensor Products**](./07-tensor-products.md) — Functions, calls, returns, and why everything gets tensored
8. [**Implementing It All**](./08-implementation.md) — Building it in Rust, the full NPA-TP pipeline

Each lecture has a problem set. Do them. You won't understand this stuff by reading.

## Philosophy

We're going to derive things rather than define them. When we need a lattice, we'll discover it by trying to solve a problem. When we need semirings, we'll notice the pattern across multiple domains. When Newton's method shows up, it'll be because we're desperate for something faster.

Math is discovered, not handed down from on high. These lectures try to recreate that experience.

## How Long Will This Take?

If you're doing the problem sets properly: maybe 20-40 hours total. You could rush through the readings in a few hours, but you wouldn't learn anything.

## Let's Go

Start with [Lecture 1: What Are We Even Doing?](./01-what-are-we-doing.md)
