use std::rc::Rc;
use crate::semiring::{Semiring, Admissible};

/// Regular expression over a semiring S
/// Uses Rc for subexpression sharing (important for efficiency)
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Expr<S: Semiring> {
    /// Additive identity (abstract zero)
    Zero,
    /// Multiplicative identity (abstract one)
    One,
    /// Constant value
    Const(S),
    /// Variable reference (index into value vector)
    Var(usize),
    /// Combine (⊕)
    Combine(Rc<Expr<S>>, Rc<Expr<S>>),
    /// Extend (⊗)
    Extend(Rc<Expr<S>>, Rc<Expr<S>>),
    /// Kleene star
    Star(Rc<Expr<S>>),
}

impl<S: Semiring> Expr<S> {
    pub fn zero() -> Self {
        Expr::Zero
    }

    pub fn one() -> Self {
        Expr::One
    }

    pub fn var(i: usize) -> Self {
        Expr::Var(i)
    }

    pub fn constant(s: S) -> Self {
        Expr::Const(s)
    }

    pub fn combine(self, other: Self) -> Self {
        // Simplify: 0 ⊕ x = x, x ⊕ 0 = x
        match (&self, &other) {
            (Expr::Zero, _) => other,
            (_, Expr::Zero) => self,
            _ => Expr::Combine(Rc::new(self), Rc::new(other)),
        }
    }

    pub fn extend(self, other: Self) -> Self {
        // Simplify: 0 ⊗ x = 0, x ⊗ 0 = 0, 1 ⊗ x = x, x ⊗ 1 = x
        match (&self, &other) {
            (Expr::Zero, _) | (_, Expr::Zero) => Expr::Zero,
            (Expr::One, _) => other,
            (_, Expr::One) => self,
            _ => Expr::Extend(Rc::new(self), Rc::new(other)),
        }
    }

    pub fn star(self) -> Self {
        // Simplify: 0* = 1, 1* = 1
        match &self {
            Expr::Zero => Expr::One,
            Expr::One => Expr::One,
            _ => Expr::Star(Rc::new(self)),
        }
    }

    /// Evaluate expression with given variable values
    pub fn eval(&self, vars: &[S]) -> S {
        match self {
            Expr::Zero => S::zero(),
            Expr::One => S::one(),
            Expr::Const(s) => s.clone(),
            Expr::Var(i) => vars[*i].clone(),
            Expr::Combine(a, b) => a.eval(vars).combine(&b.eval(vars)),
            Expr::Extend(a, b) => a.eval(vars).extend(&b.eval(vars)),
            Expr::Star(a) => a.eval(vars).star(),
        }
    }
    
    /// Evaluate with a properly-sized identity element
    /// This is needed when Expr::One appears and we need a concrete size
    pub fn eval_with_one(&self, vars: &[S], one: &S) -> S {
        match self {
            Expr::Zero => {
                // Use one's size to create a properly-sized zero
                // For now, rely on combine semantics to handle this
                S::zero()
            },
            Expr::One => one.clone(),
            Expr::Const(s) => s.clone(),
            Expr::Var(i) => vars[*i].clone(),
            Expr::Combine(a, b) => a.eval_with_one(vars, one).combine(&b.eval_with_one(vars, one)),
            Expr::Extend(a, b) => a.eval_with_one(vars, one).extend(&b.eval_with_one(vars, one)),
            Expr::Star(a) => a.eval_with_one(vars, one).star(),
        }
    }
}

/// Generalized regular expression (Definition 4.5)
/// Mixes operators from S and ST
#[derive(Clone, Debug)]
pub enum GenExpr<S: Admissible> {
    // Expressions evaluating to S
    SConst(S),
    SVar(usize),                                    // νi
    SCombine(Rc<GenExpr<S>>, Rc<GenExpr<S>>),       // ⊕
    SExtend(Rc<GenExpr<S>>, Rc<GenExpr<S>>),        // ⊗
    SStar(Rc<GenExpr<S>>),                          // *

    // Expressions evaluating to ST
    TConst(S::Tensor),
    Coupling(Rc<GenExpr<S>>, Rc<GenExpr<S>>),       // aᵗ ⊗ b (creates ST from S)
    TCombine(Rc<GenExpr<S>>, Rc<GenExpr<S>>),       // ⊕T
    TExtend(Rc<GenExpr<S>>, Rc<GenExpr<S>>),        // ⊗T
    TStar(Rc<GenExpr<S>>),                          // *T
}

/// Result of evaluating a GenExpr - either S or ST
#[derive(Clone, Debug)]
pub enum GenValue<S: Admissible> {
    S(S),
    T(S::Tensor),
}

impl<S: Admissible> GenExpr<S> {
    // S-valued constructors
    pub fn s_const(s: S) -> Self {
        GenExpr::SConst(s)
    }

    pub fn s_var(i: usize) -> Self {
        GenExpr::SVar(i)
    }

    pub fn s_combine(a: Self, b: Self) -> Self {
        GenExpr::SCombine(Rc::new(a), Rc::new(b))
    }

    pub fn s_extend(a: Self, b: Self) -> Self {
        GenExpr::SExtend(Rc::new(a), Rc::new(b))
    }

    pub fn s_star(a: Self) -> Self {
        GenExpr::SStar(Rc::new(a))
    }

    // ST-valued constructors
    pub fn t_const(t: S::Tensor) -> Self {
        GenExpr::TConst(t)
    }

    /// Coupling: C(a, b) = aᵗ ⊗ b
    pub fn coupling(a: Self, b: Self) -> Self {
        GenExpr::Coupling(Rc::new(a), Rc::new(b))
    }

    pub fn t_combine(a: Self, b: Self) -> Self {
        GenExpr::TCombine(Rc::new(a), Rc::new(b))
    }

    pub fn t_extend(a: Self, b: Self) -> Self {
        GenExpr::TExtend(Rc::new(a), Rc::new(b))
    }

    pub fn t_star(a: Self) -> Self {
        GenExpr::TStar(Rc::new(a))
    }

    /// Evaluate and return GenValue
    pub fn eval(&self, vars: &[S]) -> GenValue<S> {
        match self {
            // S-valued
            GenExpr::SConst(s) => GenValue::S(s.clone()),
            GenExpr::SVar(i) => GenValue::S(vars[*i].clone()),
            GenExpr::SCombine(a, b) => {
                let a_val = a.eval_s(vars);
                let b_val = b.eval_s(vars);
                GenValue::S(a_val.combine(&b_val))
            }
            GenExpr::SExtend(a, b) => {
                let a_val = a.eval_s(vars);
                let b_val = b.eval_s(vars);
                GenValue::S(a_val.extend(&b_val))
            }
            GenExpr::SStar(a) => {
                let a_val = a.eval_s(vars);
                GenValue::S(a_val.star())
            }

            // ST-valued
            GenExpr::TConst(t) => GenValue::T(t.clone()),
            GenExpr::Coupling(a, b) => {
                let a_val = a.eval_s(vars);
                let b_val = b.eval_s(vars);
                GenValue::T(a_val.transpose().tensor(&b_val))
            }
            GenExpr::TCombine(a, b) => {
                let a_val = a.eval_t(vars);
                let b_val = b.eval_t(vars);
                GenValue::T(a_val.combine(&b_val))
            }
            GenExpr::TExtend(a, b) => {
                let a_val = a.eval_t(vars);
                let b_val = b.eval_t(vars);
                GenValue::T(a_val.extend(&b_val))
            }
            GenExpr::TStar(a) => {
                let a_val = a.eval_t(vars);
                GenValue::T(a_val.star())
            }
        }
    }

    /// Evaluate expecting S result
    pub fn eval_s(&self, vars: &[S]) -> S {
        match self.eval(vars) {
            GenValue::S(s) => s,
            GenValue::T(_) => panic!("Expected S value, got T"),
        }
    }

    /// Evaluate expecting ST result
    pub fn eval_t(&self, vars: &[S]) -> S::Tensor {
        match self.eval(vars) {
            GenValue::T(t) => t,
            GenValue::S(_) => panic!("Expected T value, got S"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean_matrix::BoolMatrix;

    #[test]
    fn test_expr_eval() {
        // Create a simple 2x2 identity
        let id = BoolMatrix::identity(2);
        
        // Expression: var(0) ⊕ const(id)
        let expr = Expr::var(0).combine(Expr::constant(id.clone()));
        
        // With var(0) = zero, should get id
        let zero = BoolMatrix::new(2);
        let result = expr.eval(&[zero]);
        assert_eq!(result, id);
    }

    #[test]
    fn test_expr_extend() {
        let id = BoolMatrix::identity(2);
        
        // Expression: var(0) ⊗ var(0) = id ⊗ id = id
        let expr = Expr::var(0).extend(Expr::var(0));
        let result = expr.eval(&[id.clone()]);
        assert_eq!(result, id);
    }

    #[test]
    fn test_zero_one_simplification() {
        let id = BoolMatrix::identity(2);
        
        // 0 ⊕ x = x
        let expr: Expr<BoolMatrix> = Expr::zero().combine(Expr::constant(id.clone()));
        assert!(matches!(expr, Expr::Const(_)));
        
        // x ⊗ 0 = 0
        let expr: Expr<BoolMatrix> = Expr::constant(id.clone()).extend(Expr::zero());
        assert!(matches!(expr, Expr::Zero));
        
        // 1 ⊗ x = x
        let expr: Expr<BoolMatrix> = Expr::one().extend(Expr::constant(id.clone()));
        assert!(matches!(expr, Expr::Const(_)));
        
        // 0* = 1
        let expr: Expr<BoolMatrix> = Expr::zero().star();
        assert!(matches!(expr, Expr::One));
    }
}
