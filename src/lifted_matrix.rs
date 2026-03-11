use crate::semiring::{Semiring, Admissible};
use crate::sparse_matrix::{SparseBoolMatrix, SparseTensorMatrix};

/// Lifted Boolean matrix with explicit Zero and One
/// This avoids the size-0 sentinel problem by making identity elements explicit
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LiftedMatrix {
    /// Additive identity: a ⊕ Zero = a, a ⊗ Zero = Zero
    Zero,
    /// Multiplicative identity: a ⊗ One = a, a ⊕ One = a ⊕ I_n
    One,
    /// Actual sized matrix
    Sized(SparseBoolMatrix),
}

impl LiftedMatrix {
    pub fn sized(m: SparseBoolMatrix) -> Self {
        LiftedMatrix::Sized(m)
    }

    pub fn identity(n: usize) -> Self {
        LiftedMatrix::Sized(SparseBoolMatrix::identity(n))
    }

    pub fn empty(n: usize) -> Self {
        LiftedMatrix::Sized(SparseBoolMatrix::new(n))
    }
}

impl Semiring for LiftedMatrix {
    fn zero() -> Self {
        LiftedMatrix::Zero
    }

    fn one() -> Self {
        LiftedMatrix::One
    }

    fn combine(&self, other: &Self) -> Self {
        use LiftedMatrix::*;
        match (self, other) {
            // Zero is identity for combine: a ⊕ 0 = a
            (Zero, x) | (x, Zero) => x.clone(),
            
            // One ⊕ One = One (idempotent identity)
            (One, One) => One,
            
            // One ⊕ Sized(m) = identity(n) ⊕ m
            (One, Sized(m)) | (Sized(m), One) => {
                let id = SparseBoolMatrix::identity(m.n);
                Sized(id.combine(m))
            }
            
            // Sized ⊕ Sized
            (Sized(a), Sized(b)) => {
                assert_eq!(a.n, b.n);
                Sized(a.combine(b))
            }
        }
    }

    fn extend(&self, other: &Self) -> Self {
        use LiftedMatrix::*;
        match (self, other) {
            // Zero annihilates: a ⊗ 0 = 0, 0 ⊗ a = 0
            (Zero, _) | (_, Zero) => Zero,
            
            // One is identity for extend: a ⊗ 1 = a, 1 ⊗ a = a
            (One, x) | (x, One) => x.clone(),
            
            // Sized ⊗ Sized
            (Sized(a), Sized(b)) => {
                assert_eq!(a.n, b.n);
                Sized(a.extend(b))
            }
        }
    }

    fn star(&self) -> Self {
        use LiftedMatrix::*;
        match self {
            // 0* = 1 (empty sum = identity)
            Zero => One,
            
            // 1* = 1 (I* = I for idempotent)
            One => One,
            
            // m* = actual star computation
            Sized(m) => Sized(m.star()),
        }
    }
}

/// Lifted tensor matrix
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LiftedTensor {
    Zero,
    One,
    Sized(SparseTensorMatrix),
}

impl Semiring for LiftedTensor {
    fn zero() -> Self {
        LiftedTensor::Zero
    }

    fn one() -> Self {
        LiftedTensor::One
    }

    fn combine(&self, other: &Self) -> Self {
        use LiftedTensor::*;
        match (self, other) {
            (Zero, x) | (x, Zero) => x.clone(),
            (One, One) => One,
            (One, Sized(m)) | (Sized(m), One) => {
                let id = SparseTensorMatrix::identity(m.n);
                Sized(id.combine(m))
            }
            (Sized(a), Sized(b)) => {
                assert_eq!(a.n, b.n);
                Sized(a.combine(b))
            }
        }
    }

    fn extend(&self, other: &Self) -> Self {
        use LiftedTensor::*;
        match (self, other) {
            (Zero, _) | (_, Zero) => Zero,
            (One, x) | (x, One) => x.clone(),
            (Sized(a), Sized(b)) => {
                assert_eq!(a.n, b.n);
                Sized(a.extend(b))
            }
        }
    }

    fn star(&self) -> Self {
        use LiftedTensor::*;
        match self {
            Zero => One,
            One => One,
            Sized(m) => Sized(m.star()),
        }
    }
}

impl Admissible for LiftedMatrix {
    type Tensor = LiftedTensor;

    fn transpose(&self) -> Self {
        use LiftedMatrix::*;
        match self {
            Zero => Zero,
            One => One,
            Sized(m) => Sized(m.transpose()),
        }
    }

    fn tensor(&self, other: &Self) -> LiftedTensor {
        use LiftedMatrix::*;
        use LiftedTensor as T;
        match (self, other) {
            // 0 ⊗ anything = 0 in tensor land
            (Zero, _) | (_, Zero) => T::Zero,
            
            // 1 ⊗ 1 = 1
            (One, One) => T::One,
            
            // 1 ⊗ m or m ⊗ 1: need to produce identity tensor mixed with m
            // This is tricky — I ⊗ m is not just m, it's a specific tensor structure
            // For now, resolve One to actual identity matrix of m's size
            (One, Sized(m)) => {
                let id = SparseBoolMatrix::identity(m.n);
                T::Sized(id.tensor(m))
            }
            (Sized(m), One) => {
                let id = SparseBoolMatrix::identity(m.n);
                T::Sized(m.tensor(&id))
            }
            
            (Sized(a), Sized(b)) => {
                assert_eq!(a.n, b.n);
                T::Sized(a.tensor(b))
            }
        }
    }

    fn detensor_transpose(t: &LiftedTensor) -> Self {
        use LiftedTensor::*;
        use LiftedMatrix as M;
        match t {
            Zero => M::Zero,
            One => M::One,
            Sized(m) => M::Sized(SparseBoolMatrix::detensor_transpose(m)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_identity_combine() {
        let a = LiftedMatrix::identity(3);
        let zero = LiftedMatrix::Zero;
        
        assert_eq!(a.combine(&zero), a);
        assert_eq!(zero.combine(&a), a);
    }

    #[test]
    fn test_zero_annihilates_extend() {
        let a = LiftedMatrix::identity(3);
        let zero = LiftedMatrix::Zero;
        
        assert_eq!(a.extend(&zero), LiftedMatrix::Zero);
        assert_eq!(zero.extend(&a), LiftedMatrix::Zero);
    }

    #[test]
    fn test_one_identity_extend() {
        let a = LiftedMatrix::identity(3);
        let one = LiftedMatrix::One;
        
        assert_eq!(a.extend(&one), a);
        assert_eq!(one.extend(&a), a);
    }

    #[test]
    fn test_star_of_zero() {
        assert_eq!(LiftedMatrix::Zero.star(), LiftedMatrix::One);
    }

    #[test]
    fn test_star_of_one() {
        assert_eq!(LiftedMatrix::One.star(), LiftedMatrix::One);
    }

    #[test]
    fn test_semiring_laws_with_sized() {
        let a = LiftedMatrix::Sized({
            let mut m = SparseBoolMatrix::new(2);
            m.set(0, 1, true);
            m
        });
        let b = LiftedMatrix::Sized({
            let mut m = SparseBoolMatrix::new(2);
            m.set(1, 0, true);
            m
        });
        let zero = LiftedMatrix::Zero;
        let one = LiftedMatrix::One;

        // Zero annihilates
        assert_eq!(a.extend(&zero), LiftedMatrix::Zero);
        assert_eq!(zero.extend(&b), LiftedMatrix::Zero);

        // One is identity for extend
        assert_eq!(a.extend(&one), a);
        assert_eq!(one.extend(&a), a);

        // Zero is identity for combine
        assert_eq!(a.combine(&zero), a);
        assert_eq!(zero.combine(&a), a);

        // Associativity with abstract elements
        assert_eq!(
            a.extend(&one).extend(&b),
            a.extend(&one.extend(&b))
        );
    }

    #[test]
    fn test_tensor_detensor() {
        let a = LiftedMatrix::identity(2);
        let tensor = a.transpose().tensor(&a);
        let recovered = LiftedMatrix::detensor_transpose(&tensor);
        assert_eq!(recovered, a);
    }

    #[test]
    fn test_tensor_with_abstract() {
        let a = LiftedMatrix::identity(2);
        let one = LiftedMatrix::One;
        let zero = LiftedMatrix::Zero;

        // Zero ⊗ anything = Zero tensor
        assert_eq!(zero.tensor(&a), LiftedTensor::Zero);
        assert_eq!(a.tensor(&zero), LiftedTensor::Zero);

        // One ⊗ One = One tensor  
        assert_eq!(one.tensor(&one), LiftedTensor::One);
    }
}
