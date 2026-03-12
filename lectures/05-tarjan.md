# Lecture 5: Tarjan's Algorithm

We have path expressions. We know they're powerful. Now we need to compute them efficiently.

Kleene's classic algorithm constructs a regex for all pairs of nodes in O(n³). But we only care about paths from the entry node. Tarjan's algorithm exploits the structure of control flow to do better.

## The Key Insight: Dominators

A node D **dominates** node N if every path from entry to N goes through D.

**Immediate dominator** (idom): The unique closest strict dominator. Every node except entry has exactly one idom.

The immediate dominator relation forms a tree rooted at entry — the **dominator tree**.

Example:
```
CFG:                    Dominator Tree:
    1                        1
   / \                      /|\
  2   3                    2 3 4
   \ /                       |
    4                        5
    |
    5
```

Node 4 is dominated by 1 (all paths go through 1). But neither 2 nor 3 dominates 4, since you can reach 4 via either.

## Why Dominators Help

**Key lemma**: On any path to N, after the last occurrence of idom(N), all remaining nodes are dominated by idom(N).

In other words: paths "descend" the dominator tree. Once you leave a dominator, you don't come back.

This lets us decompose:

```
path(N) = path(idom(N)) ⊗ dpath(N)
```

Where `dpath(N)` is the paths from idom(N) to N through nodes dominated by idom(N).

Recursively, if D₁ → D₂ → ... → Dₖ → N is the path in the dominator tree:

```
path(N) = path(D₁) ⊗ dpath(D₂) ⊗ ... ⊗ dpath(Dₖ) ⊗ dpath(N)
```

## Computing dpath

For each node N with immediate dominator D, we need to compute `dpath(N)` = paths from D to N through D's subtree.

There are two kinds of edges to consider:

**Tree edges**: Edges from D directly to N. These contribute directly.

**Non-tree edges**: Edges from some other node M to N, where M is also dominated by D.

For non-tree edges, M must be reachable from some sibling of N in the dominator tree. We track this with a **derived graph**.

## The Derived Graph

For each idom's children (siblings in domtree), we build a small derived graph:

- Nodes: the siblings
- Edges: if there's a CFG edge from somewhere dominated by sibling S to another sibling T, add derived edge S → T

The **image** of a derived edge encodes the actual paths: from S through the domtree down to the CFG edge source, then the edge itself.

```
image(S → T) = dpath(D₂) ⊗ ... ⊗ dpath(source) ⊗ edge
```

Where D₁ → D₂ → ... → source is the domtree path from S to the edge's source.

## The Algorithm

1. **Compute dominator tree** (Lengauer-Tarjan, nearly linear)

2. **Post-order traversal** of dominator tree:
   - For each node, compute `dpath` using its children's dpaths
   - Build derived graph for siblings
   - Compute path expressions in derived graph (small, use elimination)
   - Compose tree edges, derived paths, and images

3. **Reverse post-order** to compute `path`:
   - `path(entry) = (sum of backedge paths)*`
   - `path(N) = path(idom(N)) ⊗ dpath(N)`

## Worked Example

CFG with a loop:
```
1 → 2 → 3
    ↑   ↓
    +─c─+
        ↓
        4
```
Edges: a(1→2), b(2→3), c(3→2 backedge), d(3→4)

Dominator tree:
```
1
|
2
|
3
|
4
```

Node 1 dominates everyone. Node 2 dominates 3 (because all paths to 3 go through 2). Node 3 dominates 4.

**Computing dpaths bottom-up:**

- `dpath(4) = d` (tree edge from 3)
- `dpath(3) = b` (tree edge from 2)
- `dpath(2) = a` (tree edge from 1)

**Handling the backedge:**

The backedge c (3→2) goes up the dominator tree. This contributes to a cycle: from 2, go to 3 via b, then back to 2 via c. The path `bc` can repeat any number of times.

When computing `path(2)`, we include the star of all backedge paths:
- `path(2) = path(1) ⊗ dpath(2) ⊗ (bc)*`
- `= 1 ⊗ a ⊗ (bc)* = a(bc)*`

**Final results:**
- `path(1) = 1`
- `path(2) = a(bc)*`
- `path(3) = a(bc)* ⊗ b = a(bc)*b`
- `path(4) = a(bc)*b ⊗ d = a(bc)*bd`

## Key Points

The key ideas:
1. **Decomposition**: `path(N) = path(idom(N)) ⊗ dpath(N)`
2. **Tree edges** contribute directly
3. **Non-tree edges** (including backedges) contribute through derived graph paths
4. **Stars** appear from cycles in the derived graph

The full algorithm handles all the bookkeeping. The complexity is nearly linear: O(m α(m,n)) where α is the inverse Ackermann function.

## In Code

Here's the skeleton from our implementation:

```rust
pub fn tarjan<S: Semiring>(cfg: &Cfg<S>, domtree: &DomTree) -> PathExpressions<S> {
    let mut path: HashMap<usize, Expr<S>> = HashMap::default();
    let mut dpath: HashMap<usize, Expr<S>> = HashMap::default();

    // Compute derived graph
    let (derived_edges, edge_to_derived) = compute_derived_graph(cfg, domtree);

    // Compute dpath in POST-ORDER traversal of domtree
    for node in domtree.postorder() {
        // For each child, compute dpathTree and image of derived edges
        // Then compose to get dpath[child]
    }

    // Compute path in REVERSE POST-ORDER traversal
    for node in domtree.reverse_postorder() {
        if node == cfg.entry {
            // path[entry] = (sum of backedge paths)*
        } else {
            // path[node] = path[idom(node)] ⊗ dpath[node]
        }
    }

    path
}
```

## Why This Matters

For reducible CFGs (most real programs), Tarjan's algorithm is nearly linear. This is a huge win over O(n³) Kleene.

And because path expressions encode all paths symbolically, we can evaluate them in any semiring — boolean, shortest paths, transfer functions, whatever.

The path expression is computed **once**. Then we can evaluate it with different edge values as needed.

## Connection to NPA

In NPA-TP, we use Tarjan twice:

1. **On each procedure's CFG**: Converts loops into Kleene stars, giving us equations with regular RHS

2. **On the dependence graph**: Computes how variables depend on each other across the Newton iteration

The second use is what makes the whole algorithm work — it's how we solve the linearized system efficiently.

## Problem Set 5

### Problem 5.1: Dominator Tree

Draw the dominator tree for this CFG:
```
1 → 2 → 3 → 6
↓   ↓   ↑
4 → 5 →─┘
```
(Entry is 1, edges: 1→2, 1→4, 2→3, 2→5, 4→5, 5→3, 3→6)

### Problem 5.2: Immediate Dominator

For each node in Problem 5.1, state its immediate dominator.

Which node dominates node 3? (There may be multiple dominators — list them all, then identify the immediate one.)

### Problem 5.3: Path Decomposition

Using the dominator tree from 5.1, write path(6) as a composition:
```
path(6) = path(?) ⊗ dpath(?) ⊗ ... ⊗ dpath(6)
```

### Problem 5.4: Tree vs Non-Tree Edges

For node 3 in Problem 5.1:
- Which edges to 3 are tree edges (from idom(3))?
- Which are non-tree edges?
- For non-tree edges, which sibling in the domtree dominates their source?

### Problem 5.5: Derived Graph

Build the derived graph for the children of node 1 in Problem 5.1 (i.e., the nodes whose idom is 1).

What are the derived edges?

### Problem 5.6: Loop Handling

For this CFG:
```
1 → 2 → 3
    ↑   ↓
    +───+
```

Compute path(2) step by step, showing how the backedge 3→2 contributes to the (bc)* term.

### Problem 5.7: Non-Reducible CFG

Consider this non-reducible CFG:
```
1 → 2 → 4
↓   ↑
3 →─┘
↑   ↓
+───+
```
(Edges: 1→2, 1→3, 2→4, 3→2, 2→3)

Two back edges! This is non-reducible because 2 and 3 form a strongly connected component that's not a natural loop.

How does the dominator tree look? Can you still apply Tarjan's algorithm?

### Problem 5.8: Complexity

Tarjan's algorithm is O(m α(m,n)). The derived graph computation for siblings involves O(k³) elimination where k is the number of siblings.

In the worst case, every node could be a sibling (flat dominator tree). When does this happen? What's the complexity in that case?

Why doesn't this blow up in practice for typical CFGs?

---

Next: [Lecture 6: The Interprocedural Problem](./06-interprocedural.md)
