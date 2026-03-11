use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use crate::semiring::{Semiring, Admissible};

/// Sparse Boolean matrix representation
/// Stores only the (i, j) pairs that are true
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparseBoolMatrix {
    pub n: usize,
    /// Set of (row, col) pairs where the value is true
    pub entries: HashSet<(usize, usize)>,
}

impl SparseBoolMatrix {
    pub fn new(n: usize) -> Self {
        SparseBoolMatrix {
            n,
            entries: HashSet::default(),
        }
    }

    pub fn identity(n: usize) -> Self {
        let mut m = SparseBoolMatrix::new(n);
        for i in 0..n {
            m.entries.insert((i, i));
        }
        m
    }

    pub fn get(&self, i: usize, j: usize) -> bool {
        self.entries.contains(&(i, j))
    }

    pub fn set(&mut self, i: usize, j: usize, val: bool) {
        if val {
            self.entries.insert((i, j));
        } else {
            self.entries.remove(&(i, j));
        }
    }

    /// Transpose the matrix
    pub fn transpose(&self) -> Self {
        let mut m = SparseBoolMatrix::new(self.n);
        for &(i, j) in &self.entries {
            m.entries.insert((j, i));
        }
        m
    }

    /// Number of non-zero entries
    pub fn nnz(&self) -> usize {
        self.entries.len()
    }
}

impl Semiring for SparseBoolMatrix {
    fn zero() -> Self {
        SparseBoolMatrix::new(0)
    }

    fn one() -> Self {
        SparseBoolMatrix::identity(0)
    }

    fn combine(&self, other: &Self) -> Self {
        // Handle size-0 sentinels
        if self.n == 0 { return other.clone(); }
        if other.n == 0 { return self.clone(); }
        assert_eq!(self.n, other.n);

        let mut m = SparseBoolMatrix::new(self.n);
        m.entries = self.entries.union(&other.entries).cloned().collect();
        m
    }

    fn extend(&self, other: &Self) -> Self {
        // Handle size-0 sentinels
        if self.n == 0 { return other.clone(); }
        if other.n == 0 { return self.clone(); }
        assert_eq!(self.n, other.n);

        // Build column-indexed structure for other matrix (for efficient lookup)
        let mut other_by_col: HashMap<usize, Vec<usize>> = HashMap::default();
        for &(i, j) in &other.entries {
            other_by_col.entry(i).or_default().push(j);
        }

        let mut m = SparseBoolMatrix::new(self.n);
        
        // For each (i, k) in self, find all (k, j) in other
        for &(i, k) in &self.entries {
            if let Some(cols) = other_by_col.get(&k) {
                for &j in cols {
                    m.entries.insert((i, j));
                }
            }
        }
        
        m
    }

    fn star(&self) -> Self {
        if self.n == 0 {
            return SparseBoolMatrix::identity(0);
        }

        // Transitive closure via repeated squaring with identity
        // R* = I ∪ R ∪ R² ∪ R³ ∪ ...
        // For Boolean matrices, converges when R^k = R^{k+1}
        
        let mut result = SparseBoolMatrix::identity(self.n);
        result = result.combine(self);
        
        let mut prev_size = 0;
        while result.nnz() != prev_size {
            prev_size = result.nnz();
            let squared = result.extend(&result);
            result = result.combine(&squared);
        }
        
        result
    }
}

/// Sparse tensor matrix (N² × N² but stored sparsely)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SparseTensorMatrix {
    pub n: usize,
    pub n2: usize,
    /// Set of (row, col) pairs in the N² × N² space
    pub entries: HashSet<(usize, usize)>,
}

impl SparseTensorMatrix {
    pub fn new(n: usize) -> Self {
        SparseTensorMatrix {
            n,
            n2: n * n,
            entries: HashSet::default(),
        }
    }

    pub fn identity(n: usize) -> Self {
        let mut m = SparseTensorMatrix::new(n);
        let n2 = n * n;
        for i in 0..n2 {
            m.entries.insert((i, i));
        }
        m
    }

    pub fn get(&self, i: usize, j: usize) -> bool {
        self.entries.contains(&(i, j))
    }

    pub fn set(&mut self, i: usize, j: usize, val: bool) {
        if val {
            self.entries.insert((i, j));
        } else {
            self.entries.remove(&(i, j));
        }
    }

    pub fn nnz(&self) -> usize {
        self.entries.len()
    }
}

impl Semiring for SparseTensorMatrix {
    fn zero() -> Self {
        SparseTensorMatrix::new(0)
    }

    fn one() -> Self {
        SparseTensorMatrix::identity(0)
    }

    fn combine(&self, other: &Self) -> Self {
        if self.n == 0 { return other.clone(); }
        if other.n == 0 { return self.clone(); }
        assert_eq!(self.n, other.n);

        let mut m = SparseTensorMatrix::new(self.n);
        m.entries = self.entries.union(&other.entries).cloned().collect();
        m
    }

    fn extend(&self, other: &Self) -> Self {
        if self.n == 0 { return other.clone(); }
        if other.n == 0 { return self.clone(); }
        assert_eq!(self.n, other.n);

        let mut other_by_col: HashMap<usize, Vec<usize>> = HashMap::default();
        for &(i, j) in &other.entries {
            other_by_col.entry(i).or_default().push(j);
        }

        let mut m = SparseTensorMatrix::new(self.n);
        
        for &(i, k) in &self.entries {
            if let Some(cols) = other_by_col.get(&k) {
                for &j in cols {
                    m.entries.insert((i, j));
                }
            }
        }
        
        m
    }

    fn star(&self) -> Self {
        if self.n == 0 {
            return SparseTensorMatrix::identity(0);
        }

        let mut result = SparseTensorMatrix::identity(self.n);
        result = result.combine(self);
        
        let mut prev_size = 0;
        while result.nnz() != prev_size {
            prev_size = result.nnz();
            let squared = result.extend(&result);
            result = result.combine(&squared);
        }
        
        result
    }
}

impl Admissible for SparseBoolMatrix {
    type Tensor = SparseTensorMatrix;

    fn transpose(&self) -> Self {
        SparseBoolMatrix::transpose(self)
    }

    fn tensor(&self, other: &Self) -> SparseTensorMatrix {
        if self.n == 0 { return SparseTensorMatrix::new(other.n); }
        if other.n == 0 { return SparseTensorMatrix::new(self.n); }
        assert_eq!(self.n, other.n);

        let n = self.n;
        let mut t = SparseTensorMatrix::new(n);

        // Kronecker product: only set entries where both inputs are set
        for &(a, a_prime) in &self.entries {
            for &(b, b_prime) in &other.entries {
                let i = a * n + b;
                let j = a_prime * n + b_prime;
                t.entries.insert((i, j));
            }
        }

        t
    }

    fn detensor_transpose(t: &SparseTensorMatrix) -> Self {
        let n = t.n;
        if n == 0 {
            return SparseBoolMatrix::new(0);
        }

        let mut m = SparseBoolMatrix::new(n);

        // For each entry (i, j) in tensor matrix
        // Decode: i = a'*n + b, j = a*n + b'
        // Include in result(a, b') if a' = b (diagonal constraint)
        for &(i, j) in &t.entries {
            let a_prime = i / n;
            let b = i % n;
            let a = j / n;
            let b_prime = j % n;

            // Diagonal constraint: a' = b
            if a_prime == b {
                m.entries.insert((a, b_prime));
            }
        }

        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sparse_identity() {
        let id = SparseBoolMatrix::identity(4);
        assert_eq!(id.nnz(), 4);
        assert!(id.get(0, 0));
        assert!(id.get(1, 1));
        assert!(!id.get(0, 1));
    }

    #[test]
    fn test_sparse_combine() {
        let mut a = SparseBoolMatrix::new(3);
        a.set(0, 1, true);
        
        let mut b = SparseBoolMatrix::new(3);
        b.set(1, 2, true);
        
        let c = a.combine(&b);
        assert_eq!(c.nnz(), 2);
        assert!(c.get(0, 1));
        assert!(c.get(1, 2));
    }

    #[test]
    fn test_sparse_extend() {
        let mut a = SparseBoolMatrix::new(3);
        a.set(0, 1, true);
        
        let mut b = SparseBoolMatrix::new(3);
        b.set(1, 2, true);
        
        let c = a.extend(&b);
        assert_eq!(c.nnz(), 1);
        assert!(c.get(0, 2));
    }

    #[test]
    fn test_sparse_star() {
        // Chain: 0 → 1 → 2
        let mut a = SparseBoolMatrix::new(3);
        a.set(0, 1, true);
        a.set(1, 2, true);
        
        let star = a.star();
        // Should have: identity + a + a²
        // (0,0), (1,1), (2,2) from identity
        // (0,1), (1,2) from a
        // (0,2) from a²
        assert!(star.get(0, 0));
        assert!(star.get(1, 1));
        assert!(star.get(2, 2));
        assert!(star.get(0, 1));
        assert!(star.get(1, 2));
        assert!(star.get(0, 2));
    }

    #[test]
    fn test_sparse_tensor_detensor() {
        let id = SparseBoolMatrix::identity(2);
        let tensor = id.transpose().tensor(&id);
        let recovered = SparseBoolMatrix::detensor_transpose(&tensor);
        assert_eq!(recovered, id);
    }

    #[test]
    fn test_sparse_vs_dense_equivalence() {
        use crate::boolean_matrix::BoolMatrix;
        
        // Create equivalent sparse and dense matrices
        let mut sparse = SparseBoolMatrix::new(3);
        sparse.set(0, 1, true);
        sparse.set(1, 2, true);
        
        let mut dense = BoolMatrix::new(3);
        dense.set(0, 1, true);
        dense.set(1, 2, true);
        
        // Test extend
        let sparse_ext = sparse.extend(&sparse);
        let dense_ext = dense.extend(&dense);
        
        // Convert sparse to dense for comparison
        let sparse_as_dense = sparse_to_dense(&sparse_ext);
        assert_eq!(sparse_as_dense, dense_ext);
        
        // Test star
        let sparse_star = sparse.star();
        let dense_star = dense.star();
        
        let sparse_star_as_dense = sparse_to_dense(&sparse_star);
        assert_eq!(sparse_star_as_dense, dense_star);
    }
    
    fn sparse_to_dense(s: &SparseBoolMatrix) -> crate::boolean_matrix::BoolMatrix {
        let mut d = crate::boolean_matrix::BoolMatrix::new(s.n);
        for &(i, j) in &s.entries {
            d.set(i, j, true);
        }
        d
    }
}
