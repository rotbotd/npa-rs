use crate::semiring::{Admissible, Semiring};

/// N×N Boolean matrix representing a relation over 2^|predicates| states
/// For predicates P, N = 2^|P|
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoolMatrix {
    pub n: usize,
    pub data: Vec<bool>, // row-major: data[i*n + j] = M[i,j]
}

impl BoolMatrix {
    pub fn new(n: usize) -> Self {
        BoolMatrix {
            n,
            data: vec![false; n * n],
        }
    }

    pub fn get(&self, i: usize, j: usize) -> bool {
        self.data[i * self.n + j]
    }

    pub fn set(&mut self, i: usize, j: usize, val: bool) {
        self.data[i * self.n + j] = val;
    }

    /// Identity matrix
    pub fn identity(n: usize) -> Self {
        let mut m = BoolMatrix::new(n);
        for i in 0..n {
            m.set(i, i, true);
        }
        m
    }

    /// Zero matrix of given size
    pub fn zero_n(n: usize) -> Self {
        BoolMatrix::new(n)
    }

    /// Identity matrix of given size (one)
    pub fn one_n(n: usize) -> Self {
        BoolMatrix::identity(n)
    }

    /// Matrix transpose
    pub fn transposed(&self) -> Self {
        let mut m = BoolMatrix::new(self.n);
        for i in 0..self.n {
            for j in 0..self.n {
                m.set(j, i, self.get(i, j));
            }
        }
        m
    }
}

impl Semiring for BoolMatrix {
    fn zero() -> Self {
        // Size 0 as sentinel — real usage should use zero_n/one_n
        BoolMatrix::new(0)
    }

    fn one() -> Self {
        // Size 0 as sentinel — real usage should use zero_n/one_n
        BoolMatrix::identity(0)
    }

    /// Union of relations (matrix OR)
    fn combine(&self, other: &Self) -> Self {
        assert_eq!(self.n, other.n);
        let mut m = BoolMatrix::new(self.n);
        for i in 0..self.data.len() {
            m.data[i] = self.data[i] || other.data[i];
        }
        m
    }

    /// Relational composition (matrix multiply with OR/AND)
    fn extend(&self, other: &Self) -> Self {
        assert_eq!(self.n, other.n);
        let n = self.n;
        let mut m = BoolMatrix::new(n);
        for i in 0..n {
            for j in 0..n {
                let mut val = false;
                for k in 0..n {
                    val = val || (self.get(i, k) && other.get(k, j));
                }
                m.set(i, j, val);
            }
        }
        m
    }

    /// Reflexive transitive closure
    /// a* = 1 ⊕ a ⊕ a² ⊕ ... (fixpoint)
    fn star(&self) -> Self {
        let mut result = BoolMatrix::identity(self.n);
        let mut prev = BoolMatrix::new(self.n);
        
        // Iterate until fixpoint
        while result != prev {
            prev = result.clone();
            result = result.combine(&prev.extend(self));
        }
        result
    }
}

/// N²×N² Boolean matrix (tensor product of two N×N matrices)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TensorMatrix {
    pub n: usize,       // original dimension
    pub n2: usize,      // n² (tensor dimension)
    pub data: Vec<bool>, // row-major: data[i*n2 + j]
}

impl TensorMatrix {
    pub fn new(n: usize) -> Self {
        let n2 = n * n;
        TensorMatrix {
            n,
            n2,
            data: vec![false; n2 * n2],
        }
    }

    pub fn get(&self, i: usize, j: usize) -> bool {
        self.data[i * self.n2 + j]
    }

    pub fn set(&mut self, i: usize, j: usize, val: bool) {
        self.data[i * self.n2 + j] = val;
    }

    pub fn identity(n: usize) -> Self {
        let n2 = n * n;
        let mut m = TensorMatrix::new(n);
        for i in 0..n2 {
            m.set(i, i, true);
        }
        m
    }
}

impl Semiring for TensorMatrix {
    fn zero() -> Self {
        TensorMatrix::new(0)
    }

    fn one() -> Self {
        TensorMatrix::identity(0)
    }

    fn combine(&self, other: &Self) -> Self {
        assert_eq!(self.n, other.n);
        let mut m = TensorMatrix::new(self.n);
        for i in 0..self.data.len() {
            m.data[i] = self.data[i] || other.data[i];
        }
        m
    }

    fn extend(&self, other: &Self) -> Self {
        assert_eq!(self.n, other.n);
        let n2 = self.n2;
        let mut m = TensorMatrix::new(self.n);
        for i in 0..n2 {
            for j in 0..n2 {
                let mut val = false;
                for k in 0..n2 {
                    val = val || (self.get(i, k) && other.get(k, j));
                }
                m.set(i, j, val);
            }
        }
        m
    }

    fn star(&self) -> Self {
        let mut result = TensorMatrix::identity(self.n);
        let mut prev = TensorMatrix::new(self.n);
        
        while result != prev {
            prev = result.clone();
            result = result.combine(&prev.extend(self));
        }
        result
    }
}

impl Admissible for BoolMatrix {
    type Tensor = TensorMatrix;

    fn transpose(&self) -> Self {
        self.transposed()
    }

    /// Kronecker product: (R ⊗ S)[(a-1)N + b, (a'-1)N + b'] = R(a,a') ∧ S(b,b')
    fn tensor(&self, other: &Self) -> TensorMatrix {
        assert_eq!(self.n, other.n);
        let n = self.n;
        let mut t = TensorMatrix::new(n);
        
        for a in 0..n {
            for b in 0..n {
                for a_prime in 0..n {
                    for b_prime in 0..n {
                        let i = a * n + b;
                        let j = a_prime * n + b_prime;
                        t.set(i, j, self.get(a, a_prime) && other.get(b, b_prime));
                    }
                }
            }
        }
        t
    }

    /// Detensor-transpose (equation 46):
    /// (t,·)(T(A', B, A, B')) = ∃A', B : T(A', B, A, B') ∧ A' = B
    fn detensor_transpose(t: &TensorMatrix) -> Self {
        let n = t.n;
        let mut m = BoolMatrix::new(n);
        
        // For each (a, b') in output
        for a in 0..n {
            for b_prime in 0..n {
                let mut val = false;
                // Existentially quantify over a', b with constraint a' = b
                for ab in 0..n {
                    // a' = b = ab (the diagonal constraint!)
                    let i = ab * n + ab;       // (a', b) with a' = b
                    let j = a * n + b_prime;   // (a, b')
                    val = val || t.get(i, j);
                }
                m.set(a, b_prime, val);
            }
        }
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity() {
        let id = BoolMatrix::identity(2);
        assert!(id.get(0, 0));
        assert!(!id.get(0, 1));
        assert!(!id.get(1, 0));
        assert!(id.get(1, 1));
    }

    #[test]
    fn test_transpose() {
        let mut m = BoolMatrix::new(2);
        m.set(0, 1, true);
        let t = m.transposed();
        assert!(!t.get(0, 1));
        assert!(t.get(1, 0));
    }

    #[test]
    fn test_extend() {
        // Composition of identity with itself is identity
        let id = BoolMatrix::identity(2);
        let result = id.extend(&id);
        assert_eq!(result, id);
    }

    #[test]
    fn test_tensor_detensor() {
        // (t,·)(aᵗ ⊗ b) = a ⊗ b
        let mut a = BoolMatrix::new(2);
        a.set(0, 1, true); // single edge 0 → 1
        
        let mut b = BoolMatrix::new(2);
        b.set(1, 0, true); // single edge 1 → 0
        
        let t = a.transpose().tensor(&b);
        let result = BoolMatrix::detensor_transpose(&t);
        
        // Should be a ⊗ b = composition
        let expected = a.extend(&b);
        assert_eq!(result, expected);
    }
}
