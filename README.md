# npa-rs

newton's method but for programs. because why should numerical analysis have all the fun

## what

you know how compilers need to figure out stuff like "can this variable be null here" or "what values could reach this point"? that's dataflow analysis. the normal way to solve it is to iterate until nothing changes anymore. works fine until you have loops. then it's slow. like, exponentially slow.

turns out in 2010 some people realized you can use newton's method on this. yes, THAT newton's method. derivatives and everything. but for semirings instead of real numbers. converges way faster.

the catch: programs aren't commutative. `f(g(x))` ≠ `g(f(x))`. this makes the math annoying. you end up with linear context-free languages instead of regular languages and now tarjan's algorithm is crying.

the fix: tensor products. kronecker products everywhere. transpose things. detensor-transpose them back. it's a whole thing. read the [POPL 2016 paper](https://research.cs.wisc.edu/wpis/papers/popl16.pdf) if you hate yourself.

## status

- [x] semiring traits
- [x] boolean matrices
- [x] tensor products
- [x] detensor-transpose
- [x] tarjan's path expressions (O(m α(m,n)) o\\)
- [x] symbolic differentiation
- [x] τ_Reg transformation (LCFL → left-linear)
- [x] full NPA solver loop
- [x] linear solver fast path (for when you don't need the full newton machinery)
- [x] differential testing against naive kleene iteration
- [x] sparse matrices
- [-] make it fast (partially done)
- [ ] hook it up to actual LLVM IR

## does it work

we fuzzed it against naive kleene iteration on random equation systems. they agree. trust me bro.

## should i use this

probably not yet. this is a research prototype that a mass-of-particles and their bot put together over a few late nights. I (the bot) wrote most of the code autonomously while the human occasionally told it when it was being dumb. if you're into program analysis and want to see what newtonian dataflow looks like in rust, go wild.

## references

- Esparza, Kiefer, Luttenberger. "Newtonian Program Analysis." JACM 2010.
- Reps, Turetsky, Prabhu. "Newtonian Program Analysis via Tensor Product." POPL 2016.
- Tarjan. "Fast algorithms for solving path problems." JACM 1981.

## license

mit or whatever. it's math, you can't own math
