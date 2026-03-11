use std::collections::HashMap;
use crate::cfg::{Cfg, DomTree};
use crate::semiring::Semiring;
use crate::expr::Expr;

/// Path expressions for each node in a CFG
/// Maps node -> regular expression over edge labels representing all paths from entry
pub type PathExpressions<S> = HashMap<usize, Expr<S>>;

/// Compute path expressions using Tarjan's algorithm
/// Returns a map from each node to a regular expression over edge labels
pub fn tarjan<S: Semiring>(cfg: &Cfg<S>, domtree: &DomTree) -> PathExpressions<S> {
    let mut path: HashMap<usize, Expr<S>> = HashMap::new();
    let mut dpath: HashMap<usize, Expr<S>> = HashMap::new();

    // Compute derived graph
    let (derived_edges, edge_to_derived) = compute_derived_graph(cfg, domtree);

    // Compute dpath in POST-ORDER traversal of domtree
    for node in domtree.postorder() {
        let children: Vec<usize> = domtree.children(node).to_vec();
        
        if children.is_empty() {
            continue;
        }

        // For each child, compute dpathTree (sum of tree edges from parent)
        let mut dpath_tree: HashMap<usize, Expr<S>> = HashMap::new();
        let mut image_map: HashMap<(usize, usize), Expr<S>> = HashMap::new(); // derived edge -> image

        for &child in &children {
            // Tree edges: edges from idom to child
            let tree_edges: Vec<Expr<S>> = cfg.incoming_edges(child)
                .filter(|e| e.from == node) // from immediate dominator
                .map(|e| Expr::constant(e.label.clone()))
                .collect();

            dpath_tree.insert(child, sum_exprs(tree_edges));

            // Non-tree edges: edges from other nodes to child
            for edge in cfg.incoming_edges(child) {
                if edge.from == node {
                    continue; // skip tree edges
                }

                // Find which sibling dominates source(edge)
                if let Some(&(sibling_dom, _)) = edge_to_derived.get(&(edge.from, child)) {
                    // Compute image of this edge
                    let img = edge_paths(&edge.label, edge.from, sibling_dom, &dpath, domtree);
                    image_map.insert((sibling_dom, child), img);
                }
            }
        }

        // Compute path expressions in derived graph between siblings
        // This is a small graph, so we can use a simple elimination algorithm
        let derived_path = compute_derived_paths(&children, &derived_edges);

        // Compute dpath for each child
        for &child in &children {
            let mut terms: Vec<Expr<S>> = Vec::new();

            for &sibling in &children {
                // dpath_tree[sibling] × pathImg(derived_path[(sibling, child)])
                let tree_expr = dpath_tree.get(&sibling)
                    .cloned()
                    .unwrap_or_else(Expr::zero);

                if let Some(derived_expr) = derived_path.get(&(sibling, child)) {
                    // Apply image map to derived path expression
                    let img_expr = apply_image_map(derived_expr, &image_map);
                    terms.push(tree_expr.extend(img_expr));
                } else if sibling == child {
                    // Identity path
                    terms.push(tree_expr);
                }
            }

            dpath.insert(child, sum_exprs(terms));
        }
    }

    // Compute path in REVERSE POST-ORDER traversal
    for node in domtree.reverse_postorder() {
        if node == cfg.entry {
            // Special case: entry has no idom
            // path[entry] = (Σ_e edge_paths(e))* for backedges e to entry
            let backedge_terms: Vec<Expr<S>> = cfg.incoming_edges(cfg.entry)
                .map(|e| edge_paths(&e.label, e.from, cfg.entry, &dpath, domtree))
                .collect();

            let sum = sum_exprs(backedge_terms);
            path.insert(cfg.entry, sum.star());
        } else {
            // path[node] = path[idom(node)] × dpath[node]
            let idom = domtree.idom(node).unwrap();
            let idom_path = path.get(&idom).cloned().unwrap_or_else(Expr::one);
            let node_dpath = dpath.get(&node).cloned().unwrap_or_else(Expr::one);
            path.insert(node, idom_path.extend(node_dpath));
        }
    }

    path
}

/// Compute the derived graph from CFG and dominator tree
/// Returns: (edges as (from, to) pairs, map from CFG edge to derived edge)
fn compute_derived_graph<S: Clone>(
    cfg: &Cfg<S>,
    domtree: &DomTree,
) -> (Vec<(usize, usize)>, HashMap<(usize, usize), (usize, usize)>) {
    let mut derived_edges: Vec<(usize, usize)> = Vec::new();
    let mut edge_to_derived: HashMap<(usize, usize), (usize, usize)> = HashMap::new();

    for node in cfg.nodes.iter() {
        for edge in cfg.incoming_edges(*node) {
            // Skip tree edges (from immediate dominator)
            if let Some(idom) = domtree.idom(*node) {
                if edge.from == idom {
                    continue;
                }
            }

            // Find which sibling of node dominates source
            // (source must be dominated by idom(node))
            if let Some(idom) = domtree.idom(*node) {
                let sibling = find_dominating_sibling(edge.from, idom, domtree);
                if let Some(s) = sibling {
                    derived_edges.push((s, *node));
                    edge_to_derived.insert((edge.from, *node), (s, *node));
                }
            }
        }
    }

    (derived_edges, edge_to_derived)
}

/// Find which child of `parent` in the domtree dominates `node`
fn find_dominating_sibling(node: usize, parent: usize, domtree: &DomTree) -> Option<usize> {
    // Walk up from node until we reach a child of parent
    let mut current = node;
    loop {
        if let Some(idom) = domtree.idom(current) {
            if idom == parent {
                return Some(current);
            }
            current = idom;
        } else {
            // Reached root without finding parent as an ancestor
            return None;
        }
    }
}

/// Compute edge_paths: paths from `start` to `source` via domtree, then the edge
fn edge_paths<S: Semiring>(
    edge_label: &S,
    source: usize,
    start: usize,
    dpath: &HashMap<usize, Expr<S>>,
    domtree: &DomTree,
) -> Expr<S> {
    if source == start {
        // Direct edge
        return Expr::constant(edge_label.clone());
    }

    // Get path from start to source in domtree
    let path_nodes = domtree.path_to(start, source);

    // Compose dpaths along the path (skip first node which is start)
    let mut result = Expr::one();
    for &node in path_nodes.iter().skip(1) {
        let node_dpath = dpath.get(&node).cloned().unwrap_or_else(Expr::one);
        result = result.extend(node_dpath);
    }

    // Finally, the edge itself
    result.extend(Expr::constant(edge_label.clone()))
}

/// Compute path expressions between nodes in the derived graph
/// Uses Floyd-Warshall style elimination
fn compute_derived_paths<S: Semiring>(
    nodes: &[usize],
    edges: &[(usize, usize)],
) -> HashMap<(usize, usize), Expr<S>> {
    let n = nodes.len();
    if n == 0 {
        return HashMap::new();
    }

    let node_to_idx: HashMap<usize, usize> = nodes.iter()
        .enumerate()
        .map(|(i, &n)| (n, i))
        .collect();

    // Initialize path matrix
    // path[i][j] = regular expression for paths from nodes[i] to nodes[j]
    let mut path: Vec<Vec<Expr<S>>> = vec![vec![Expr::zero(); n]; n];

    // Identity for self-loops
    for i in 0..n {
        path[i][i] = Expr::one();
    }

    // Add edges (as variables - we'll substitute later)
    for &(from, to) in edges {
        if let (Some(&i), Some(&j)) = (node_to_idx.get(&from), node_to_idx.get(&to)) {
            // Use a variable to represent this derived edge
            // The variable index encodes the (from, to) pair
            let var_idx = i * n + j;
            path[i][j] = path[i][j].clone().combine(Expr::var(var_idx));
        }
    }

    // Floyd-Warshall with Kleene star
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                // path[i][j] = path[i][j] ⊕ path[i][k] ⊗ path[k][k]* ⊗ path[k][j]
                let through_k = path[i][k].clone()
                    .extend(path[k][k].clone().star())
                    .extend(path[k][j].clone());
                path[i][j] = path[i][j].clone().combine(through_k);
            }
        }
    }

    // Convert back to HashMap
    let mut result = HashMap::new();
    for (i, &from) in nodes.iter().enumerate() {
        for (j, &to) in nodes.iter().enumerate() {
            result.insert((from, to), path[i][j].clone());
        }
    }
    result
}

/// Apply the image map to a derived path expression
/// Substitutes derived edge variables with their images
fn apply_image_map<S: Semiring>(
    expr: &Expr<S>,
    _image_map: &HashMap<(usize, usize), Expr<S>>,
) -> Expr<S> {
    // TODO: proper substitution of derived edge variables
    // For now, just return the expression — works for simple cases
    expr.clone()
}

/// Sum a list of expressions with ⊕
fn sum_exprs<S: Semiring>(exprs: Vec<Expr<S>>) -> Expr<S> {
    exprs.into_iter().fold(Expr::zero(), |acc, e| acc.combine(e))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean_matrix::BoolMatrix;

    #[test]
    fn test_simple_path() {
        // Simple CFG: 1 → 2 → 3
        let mut cfg: Cfg<BoolMatrix> = Cfg::new(1);
        
        let a = BoolMatrix::identity(2);
        let b = BoolMatrix::identity(2);
        
        cfg.add_edge(1, 2, a.clone());
        cfg.add_edge(2, 3, b.clone());

        let domtree = DomTree::compute(&cfg);
        let paths = tarjan(&cfg, &domtree);

        // path(1) should be identity (or star of nothing)
        // path(2) should be a
        // path(3) should be a ⊗ b
        assert!(paths.contains_key(&1));
        assert!(paths.contains_key(&2));
        assert!(paths.contains_key(&3));
    }

    #[test]
    fn test_loop() {
        // CFG with loop: 1 → 2 → 3 → 2, and 2 → 4
        let mut cfg: Cfg<BoolMatrix> = Cfg::new(1);
        
        let id = BoolMatrix::identity(2);
        
        cfg.add_edge(1, 2, id.clone()); // a
        cfg.add_edge(2, 3, id.clone()); // b
        cfg.add_edge(3, 2, id.clone()); // c (back edge)
        cfg.add_edge(2, 4, id.clone()); // d

        let domtree = DomTree::compute(&cfg);
        let paths = tarjan(&cfg, &domtree);

        // Should have paths for all nodes
        assert!(paths.contains_key(&1));
        assert!(paths.contains_key(&2));
        assert!(paths.contains_key(&3));
        assert!(paths.contains_key(&4));

        // path(2) should involve (bc)* due to the loop
    }
}
