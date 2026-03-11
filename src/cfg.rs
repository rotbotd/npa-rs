use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

/// A labeled edge in the CFG
#[derive(Clone, Debug)]
pub struct Edge<L> {
    pub from: usize,
    pub to: usize,
    pub label: L,
}

/// Control Flow Graph with labeled edges
#[derive(Clone, Debug)]
pub struct Cfg<L> {
    pub entry: usize,
    pub nodes: HashSet<usize>,
    pub edges: Vec<Edge<L>>,
    // Adjacency lists for efficiency
    pub succs: HashMap<usize, Vec<usize>>,  // node -> successor nodes
    pub preds: HashMap<usize, Vec<usize>>,  // node -> predecessor nodes
    pub out_edges: HashMap<usize, Vec<usize>>,  // node -> edge indices
    pub in_edges: HashMap<usize, Vec<usize>>,   // node -> edge indices
}

impl<L: Clone> Cfg<L> {
    pub fn new(entry: usize) -> Self {
        let mut nodes = HashSet::default();
        nodes.insert(entry);
        Cfg {
            entry,
            nodes,
            edges: Vec::new(),
            succs: HashMap::default(),
            preds: HashMap::default(),
            out_edges: HashMap::default(),
            in_edges: HashMap::default(),
        }
    }

    pub fn add_node(&mut self, node: usize) {
        self.nodes.insert(node);
    }

    pub fn add_edge(&mut self, from: usize, to: usize, label: L) {
        self.nodes.insert(from);
        self.nodes.insert(to);
        
        let edge_idx = self.edges.len();
        self.edges.push(Edge { from, to, label });
        
        self.succs.entry(from).or_default().push(to);
        self.preds.entry(to).or_default().push(from);
        self.out_edges.entry(from).or_default().push(edge_idx);
        self.in_edges.entry(to).or_default().push(edge_idx);
    }

    pub fn successors(&self, node: usize) -> &[usize] {
        self.succs.get(&node).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn predecessors(&self, node: usize) -> &[usize] {
        self.preds.get(&node).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn incoming_edges(&self, node: usize) -> impl Iterator<Item = &Edge<L>> {
        self.in_edges
            .get(&node)
            .map(|indices| indices.as_slice())
            .unwrap_or(&[])
            .iter()
            .map(|&i| &self.edges[i])
    }

    pub fn outgoing_edges(&self, node: usize) -> impl Iterator<Item = &Edge<L>> {
        self.out_edges
            .get(&node)
            .map(|indices| indices.as_slice())
            .unwrap_or(&[])
            .iter()
            .map(|&i| &self.edges[i])
    }
}

/// Dominator tree
#[derive(Clone, Debug)]
pub struct DomTree {
    pub entry: usize,
    pub idom: HashMap<usize, usize>,           // node -> immediate dominator
    pub children: HashMap<usize, Vec<usize>>,  // node -> dominated children
}

impl DomTree {
    /// Compute dominator tree using simple iterative algorithm
    /// (Not Lengauer-Tarjan, but correct and simple)
    pub fn compute<L: Clone>(cfg: &Cfg<L>) -> Self {
        let nodes: Vec<usize> = cfg.nodes.iter().copied().collect();
        let n = nodes.len();
        
        // Map nodes to indices for the algorithm
        let node_to_idx: HashMap<usize, usize> = nodes.iter()
            .enumerate()
            .map(|(i, &n)| (n, i))
            .collect();
        
        // Initialize dominators: entry dominates itself, others dominated by all
        let mut doms: Vec<HashSet<usize>> = vec![HashSet::default(); n];
        let entry_idx = node_to_idx[&cfg.entry];
        doms[entry_idx].insert(cfg.entry);
        
        for (i, &node) in nodes.iter().enumerate() {
            if node != cfg.entry {
                doms[i] = cfg.nodes.clone();
            }
        }
        
        // Iterate until fixpoint
        let mut changed = true;
        while changed {
            changed = false;
            for &node in &nodes {
                if node == cfg.entry {
                    continue;
                }
                let idx = node_to_idx[&node];
                
                // Dom(n) = {n} ∪ ⋂_{p ∈ preds(n)} Dom(p)
                let preds = cfg.predecessors(node);
                if preds.is_empty() {
                    continue;
                }
                
                let mut new_dom: HashSet<usize> = doms[node_to_idx[&preds[0]]].clone();
                for &pred in &preds[1..] {
                    new_dom = new_dom.intersection(&doms[node_to_idx[&pred]])
                        .copied()
                        .collect();
                }
                new_dom.insert(node);
                
                if new_dom != doms[idx] {
                    doms[idx] = new_dom;
                    changed = true;
                }
            }
        }
        
        // Extract immediate dominators
        let mut idom = HashMap::default();
        for &node in &nodes {
            if node == cfg.entry {
                continue;
            }
            let idx = node_to_idx[&node];
            
            // idom(n) = the dominator of n that is dominated by all other dominators of n
            // (i.e., the "closest" strict dominator)
            let strict_doms: HashSet<usize> = doms[idx].iter()
                .filter(|&&d| d != node)
                .copied()
                .collect();
            
            for &candidate in &strict_doms {
                let candidate_idx = node_to_idx[&candidate];
                // candidate is idom if it dominates all other strict dominators
                let candidate_doms = &doms[candidate_idx];
                let is_idom = strict_doms.iter()
                    .all(|&d| d == candidate || candidate_doms.contains(&d));
                if is_idom {
                    idom.insert(node, candidate);
                    break;
                }
            }
        }
        
        // Build children map
        let mut children: HashMap<usize, Vec<usize>> = HashMap::default();
        for &node in &nodes {
            children.insert(node, Vec::new());
        }
        for (&node, &parent) in &idom {
            children.get_mut(&parent).unwrap().push(node);
        }
        
        DomTree {
            entry: cfg.entry,
            idom,
            children,
        }
    }

    pub fn idom(&self, node: usize) -> Option<usize> {
        self.idom.get(&node).copied()
    }

    pub fn children(&self, node: usize) -> &[usize] {
        self.children.get(&node).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Post-order traversal of dominator tree
    pub fn postorder(&self) -> Vec<usize> {
        let mut result = Vec::new();
        let mut stack = vec![(self.entry, false)];
        
        while let Some((node, visited)) = stack.pop() {
            if visited {
                result.push(node);
            } else {
                stack.push((node, true));
                for &child in self.children(node) {
                    stack.push((child, false));
                }
            }
        }
        result
    }

    /// Reverse post-order (topological order from entry)
    pub fn reverse_postorder(&self) -> Vec<usize> {
        let mut result = self.postorder();
        result.reverse();
        result
    }

    /// Get path from ancestor to descendant in domtree
    pub fn path_to(&self, from: usize, to: usize) -> Vec<usize> {
        // Walk up from `to` until we reach `from`
        let mut path = Vec::new();
        let mut current = to;
        while current != from {
            path.push(current);
            current = self.idom(current).expect("to must be dominated by from");
        }
        path.push(from);
        path.reverse();
        path
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cfg_basic() {
        let mut cfg = Cfg::new(1);
        cfg.add_edge(1, 2, "a");
        cfg.add_edge(2, 3, "b");
        cfg.add_edge(3, 2, "c");  // back edge
        cfg.add_edge(2, 4, "d");
        
        assert_eq!(cfg.nodes.len(), 4);
        assert_eq!(cfg.successors(2), &[3, 4]);
        assert_eq!(cfg.predecessors(2), &[1, 3]);
    }

    #[test]
    fn test_domtree() {
        let mut cfg = Cfg::new(1);
        cfg.add_edge(1, 2, "a");
        cfg.add_edge(2, 3, "b");
        cfg.add_edge(3, 2, "c");  // back edge
        cfg.add_edge(2, 4, "d");
        
        let domtree = DomTree::compute(&cfg);
        
        // 1 dominates everything
        // 2 dominates 3 and 4
        assert_eq!(domtree.idom(2), Some(1));
        assert_eq!(domtree.idom(3), Some(2));
        assert_eq!(domtree.idom(4), Some(2));
        
        assert!(domtree.children(1).contains(&2));
        assert!(domtree.children(2).contains(&3));
        assert!(domtree.children(2).contains(&4));
    }

    #[test]
    fn test_domtree_postorder() {
        let mut cfg = Cfg::new(1);
        cfg.add_edge(1, 2, "a");
        cfg.add_edge(2, 3, "b");
        cfg.add_edge(2, 4, "d");
        
        let domtree = DomTree::compute(&cfg);
        let postorder = domtree.postorder();
        
        // Children before parents
        let pos = |n: usize| postorder.iter().position(|&x| x == n).unwrap();
        assert!(pos(3) < pos(2));
        assert!(pos(4) < pos(2));
        assert!(pos(2) < pos(1));
    }
}
