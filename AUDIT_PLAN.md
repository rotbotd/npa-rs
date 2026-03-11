# NPA-TP Audit Plan

Tracking all issues raised in Gemini's audit. Items marked with ✅ are resolved.

## Phase 1: `semiring.rs` and `expr.rs`

### 1A. Compile-time Type Safety for GenExpr
**Status:** 🔲 Not started

Split `GenExpr` into two mutually recursive types (`GenExpr` for S, `GenExprT` for S::Tensor) to eliminate runtime panics from `eval_s()`/`eval_t()` type mismatches.

**Priority:** Low - current code works, this is a refinement

### 1B. Subexpression Sharing / Hash Consing
**Status:** 🔲 Not started

Current `Rc` usage doesn't actually deduplicate identical subtrees. Consider hash consing or arena allocation (`bumpalo`) if memory becomes an issue during Tarjan's algorithm.

**Priority:** Low - optimize only if needed

### 1C. Constant Folding in GenExpr
**Status:** 🔲 Not started

`Expr` has constant folding (`0 ⊕ x = x`, `0* = 1`) but `GenExpr` doesn't. Should implement the same simplifications.

**Priority:** Low - refinement

---

## Phase 2: `boolean_matrix.rs`, `sparse_matrix.rs`, `lifted_matrix.rs`

### 2A. HashSet Overhead in sparse_matrix.rs
**Status:** 🔲 Not started

Replace `HashSet`/`HashMap` with `FxHashSet`/`FxHashMap` from `rustc-hash` crate for ~2x speedup on integer keys.

**Priority:** Medium - easy win

### 2B. Memory Allocations in sparse_matrix::extend
**Status:** 🔲 Not started

Currently allocates new `HashMap` and `Vec`s every iteration. Consider switching to `Vec<Vec<usize>>` adjacency list representation.

**Priority:** Medium - performance

### 2C. Dense Matrix Bit Packing
**Status:** 🔲 Not started

`Vec<bool>` uses 1 byte per boolean. Use `bitvec` or `fixedbitset` for 8x memory reduction and SIMD bitwise operations.

**Priority:** Low - optimize if needed

### 2D. Repeated Squaring for Dense star()
**Status:** 🔲 Not started

`BoolMatrix::star` uses linear iteration O(k), but `SparseBoolMatrix::star` uses repeated squaring O(log k). Update dense to match.

**Priority:** Low - performance

---

## Phase 3: `cfg.rs` and `tarjan.rs`

### 3A. Overwriting image_map Bug
**Status:** ✅ Resolved

Multiple non-tree edges mapping to the same derived edge would overwrite each other. Fixed by using `entry().and_modify().or_insert()` to combine images.

**Fix:** Commit `8e132b1` - "fix: thread size through NPA to fix sentinel zero bug" (included Tarjan fix)

### 3B. Derived Edges Deduplication
**Status:** 🔲 Not started

`derived_edges` Vec can contain duplicates, causing redundant `X ⊕ X` in AST. Should deduplicate before Floyd-Warshall.

**Priority:** Low - mathematically correct, just inefficient

### 3C. FxHash for cfg.rs and tarjan.rs
**Status:** 🔲 Not started

Same as 2A - replace standard HashMap/HashSet with FxHash variants.

**Priority:** Medium - easy win

---

## Phase 4: `differentiate.rs` and `regularize.rs`

### 4A. X vs Y Shadowing Bug - Direct Coefficient Tracking
**Status:** ✅ Resolved

The old `differentiate` + `extract_lcfl_terms` approach confused the linear placeholder Y with the constant X evaluated at ν. Implemented `differentiate_to_lcfl` which computes LCFL coefficients directly during recursive descent.

**Fix:** Implemented in `regularize.rs` - handles Star without panicking

### 4B. Delete Old differentiate/extract Code
**Status:** 🔲 Not started

Now that `differentiate_to_lcfl` exists, the old `differentiate` function and `extract_lcfl_terms` can be removed or deprecated.

**Priority:** Low - cleanup

---

## Phase 5: `npa.rs`

### 5A. Missing Zero/One Variants in eval_path_expr
**Status:** ✅ Resolved

Pattern match wasn't exhaustive - added `Expr::Zero` and `Expr::One` cases with proper sized identity/zero handling.

**Fix:** Commit `8e132b1`

### 5B. Hoist Formal Derivative Out of Loop
**Status:** 🔲 Not started

Currently computes `differentiate_to_lcfl` inside the Newton loop. The derivative AST never changes between rounds - only the evaluation at ν changes. Should precompute in `NpaSolver::new()`.

**Priority:** Medium - performance

### 5C. Memoization for eval_path_expr
**Status:** 🔲 Not started

Path expressions form a DAG but `eval_path_expr` traverses as if it were a tree, causing O(2^N) blowup. Use `Rc` pointer address as memoization key.

**Priority:** High - critical for large CFGs

### 5D. Size Threading for Sentinel Zeros
**Status:** ✅ Resolved

`S::zero()` returned size-0 sentinel matrices which were incorrectly treated as identity in `extend`. Fixed by:
- Adding `zero_sized(n)`, `one_sized(n)`, `size()` to Semiring trait
- NPA uses properly-sized zeros for initialization
- `eval_path_expr` takes `tensor_size` parameter

**Fix:** Commit `8e132b1`

---

## Phase 6: `fuzz.rs`

### 6A. Linearity Trap - Use Non-Linear Expressions
**Status:** ✅ Resolved

`test_npa_vs_naive` was using `random_linear_expr` which only tests trivial cases. Changed to use `random_expr` to stress test quadratic equations and stars.

**Fix:** Commit `1c7aba5` - "fuzz: use non-linear expressions and add detensor-transpose law"

### 6B. Detensor-Transpose Law Test
**Status:** ✅ Resolved

Added the fundamental theorem (Equation 43) to `test_admissible_laws`:
```
(t,·)(aᵀ ⊗ b) = a ⊗ b
```

**Fix:** Commit `1c7aba5`

---

## Summary

| Phase | Total Issues | Resolved | Remaining |
|-------|-------------|----------|-----------|
| 1     | 3           | 0        | 3         |
| 2     | 4           | 0        | 4         |
| 3     | 3           | 1        | 2         |
| 4     | 2           | 1        | 1         |
| 5     | 4           | 3        | 1         |
| 6     | 2           | 2        | 0         |
| **Total** | **18**  | **7**    | **11**    |

### High Priority Remaining
- **5C.** Memoization for eval_path_expr (O(2^N) → O(N))

### Medium Priority Remaining  
- **2A, 3C.** FxHash everywhere (easy 2x speedup)
- **2B.** Sparse matrix allocation in hot path
- **5B.** Hoist derivative computation

### Low Priority (Refinements)
- 1A, 1B, 1C, 2C, 2D, 3B, 4B
