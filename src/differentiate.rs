use std::rc::Rc;
use crate::semiring::Semiring;
use crate::expr::Expr;

/// Compute the differential D_Xj[expr] with respect to variable j
/// 
/// The differential is a linear operator that extracts the "coefficient" of Xj
/// in the expression, treating Xj as a formal variable.
/// 
/// Key rules:
/// - D_Xj[c] = 0                           (constants)
/// - D_Xj[Xi] = 1 if i=j, else 0           (variables)
/// - D_Xj[a ⊕ b] = D_Xj[a] ⊕ D_Xj[b]       (combine)
/// - D_Xj[a ⊗ b] = D_Xj[a] ⊗ b ⊕ a ⊗ D_Xj[b]  (product rule - BUT non-commutative!)
/// - D_Xj[a*] = a* ⊗ D_Xj[a] ⊗ a*          (Kleene star - Theorem 6.3)
/// 
/// For non-commutative semirings, the product rule gives us an LCFL, not regular.
/// That's why we need τ_Reg transformation after differentiation.
pub fn differentiate<S: Semiring>(expr: &Expr<S>, var_idx: usize) -> Expr<S> {
    match expr {
        Expr::Const(_) => Expr::zero(),
        
        Expr::Var(i) => {
            if *i == var_idx {
                Expr::one()
            } else {
                Expr::zero()
            }
        }
        
        Expr::Combine(a, b) => {
            let da = differentiate(a, var_idx);
            let db = differentiate(b, var_idx);
            da.combine(db)
        }
        
        Expr::Extend(a, b) => {
            // Product rule: D[a ⊗ b] = D[a] ⊗ b ⊕ a ⊗ D[b]
            // This creates LCFL terms like D[a] ⊗ b (left-linear) 
            // and a ⊗ D[b] (right-linear)
            let da = differentiate(a, var_idx);
            let db = differentiate(b, var_idx);
            
            // D[a] ⊗ b
            let left = Expr::Extend(Rc::new(da), b.clone());
            // a ⊗ D[b]
            let right = Expr::Extend(a.clone(), Rc::new(db));
            
            left.combine(right)
        }
        
        Expr::Star(a) => {
            // D[a*] = a* ⊗ D[a] ⊗ a* (Theorem 6.3)
            // The result is still linear in Y because D[a] is linear
            let da = differentiate(a, var_idx);
            let a_star = Expr::Star(a.clone());
            
            // a* ⊗ D[a] ⊗ a*
            a_star.clone().extend(da).extend(a_star)
        }
    }
}

/// Evaluate an expression at a given point, substituting variable values
/// This is used to get ν-specific coefficients: D_Xj[f]|ν
pub fn eval_at<S: Semiring>(expr: &Expr<S>, values: &[S]) -> S {
    expr.eval(values)
}

/// Extract the linear coefficient structure from a differentiated expression
/// 
/// After differentiation, we have terms like:
/// - c ⊗ Y ⊗ d  (Y appears once, sandwiched)
/// - c ⊗ Y      (Y on right)
/// - Y ⊗ d      (Y on left)
/// - c          (no Y - shouldn't happen if properly differentiated)
/// 
/// This is the LCFL structure that τ_Reg will transform.
#[derive(Clone, Debug)]
pub struct LinearTerm<S: Semiring> {
    pub left: Expr<S>,   // coefficient on left of Y
    pub right: Expr<S>,  // coefficient on right of Y
}

impl<S: Semiring> LinearTerm<S> {
    pub fn new(left: Expr<S>, right: Expr<S>) -> Self {
        LinearTerm { left, right }
    }
    
    /// Create term: 1 ⊗ Y ⊗ 1 (just Y)
    pub fn identity() -> Self {
        LinearTerm {
            left: Expr::one(),
            right: Expr::one(),
        }
    }
    
    /// Create term: c ⊗ Y ⊗ 1 (Y on right)
    pub fn left_coeff(c: Expr<S>) -> Self {
        LinearTerm {
            left: c,
            right: Expr::one(),
        }
    }
    
    /// Create term: 1 ⊗ Y ⊗ c (Y on left)
    pub fn right_coeff(c: Expr<S>) -> Self {
        LinearTerm {
            left: Expr::one(),
            right: c,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean_matrix::BoolMatrix;
    
    #[test]
    fn test_diff_constant() {
        let c = Expr::<BoolMatrix>::constant(BoolMatrix::identity(2));
        let dc = differentiate(&c, 0);
        
        // Derivative of constant is zero
        match dc {
            Expr::Const(m) => assert_eq!(m, BoolMatrix::new(0)),
            _ => panic!("Expected Const"),
        }
    }
    
    #[test]
    fn test_diff_var() {
        let x0 = Expr::<BoolMatrix>::var(0);
        let x1 = Expr::<BoolMatrix>::var(1);
        
        // D_X0[X0] = 1
        let dx0 = differentiate(&x0, 0);
        match dx0 {
            Expr::Const(m) => assert_eq!(m, BoolMatrix::identity(0)),
            _ => panic!("Expected one"),
        }
        
        // D_X0[X1] = 0
        let dx1 = differentiate(&x1, 0);
        match dx1 {
            Expr::Const(m) => assert_eq!(m, BoolMatrix::new(0)),
            _ => panic!("Expected zero"),
        }
    }
    
    #[test]
    fn test_diff_product() {
        // D_X0[X0 ⊗ X0] = D_X0[X0] ⊗ X0 ⊕ X0 ⊗ D_X0[X0]
        //              = 1 ⊗ X0 ⊕ X0 ⊗ 1
        //              = X0 ⊕ X0
        let x0 = Expr::<BoolMatrix>::var(0);
        let expr = x0.clone().extend(x0.clone());
        let d = differentiate(&expr, 0);
        
        // Should be Combine(Extend(one, var(0)), Extend(var(0), one))
        match d {
            Expr::Combine(_, _) => (), // Good, it's a sum
            _ => panic!("Expected Combine"),
        }
    }
    
    #[test]
    fn test_diff_star() {
        // D_X0[X0*] = X0* ⊗ D_X0[X0] ⊗ X0* = X0* ⊗ 1 ⊗ X0* = X0* ⊗ X0*
        let x0 = Expr::<BoolMatrix>::var(0);
        let expr = x0.star();
        let d = differentiate(&expr, 0);
        
        // Should be Extend(Extend(Star(var(0)), one), Star(var(0)))
        match d {
            Expr::Extend(_, _) => (), // Good
            _ => panic!("Expected Extend"),
        }
    }
}
