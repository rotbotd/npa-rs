/// A semiring with combine (⊕), extend (⊗), zero, one, and Kleene star
pub trait Semiring: Clone + Eq + std::fmt::Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn combine(&self, other: &Self) -> Self; // ⊕
    fn extend(&self, other: &Self) -> Self;  // ⊗
    fn star(&self) -> Self;                  // Kleene star: a* = 1 ⊕ a ⊕ aa ⊕ ...
    
    /// Zero with explicit size (for path expression evaluation)
    /// Default implementation just calls zero()
    fn zero_sized(_size: usize) -> Self { Self::zero() }
    
    /// One with explicit size (for path expression evaluation)  
    /// Default implementation just calls one()
    fn one_sized(_size: usize) -> Self { Self::one() }
    
    /// Get the size of this element (for matrices, returns n)
    /// Default returns 0 for unsized elements
    fn size(&self) -> usize { 0 }
}

/// An admissible semiring supports NPA-TP via tensor products
pub trait Admissible: Semiring {
    type Tensor: Semiring;

    /// Transpose: reverses order under extend
    /// (a ⊗ b)ᵗ = bᵗ ⊗ aᵗ
    fn transpose(&self) -> Self;

    /// Tensor product: S × S → ST (Kronecker product for matrices)
    fn tensor(&self, other: &Self) -> Self::Tensor;

    /// Detensor-transpose: ST → S
    /// (t,·)(aᵗ ⊗ b) = a ⊗ b
    /// Distributes over combine: (t,·)(p₁ ⊕T p₂) = (t,·)(p₁) ⊕ (t,·)(p₂)
    fn detensor_transpose(t: &Self::Tensor) -> Self;
}

/// Coupling operation: C(a, b) = aᵗ ⊗ b
pub fn coupling<S: Admissible>(a: &S, b: &S) -> S::Tensor {
    a.transpose().tensor(b)
}
