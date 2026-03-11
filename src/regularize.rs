use crate::semiring::{Semiring, Admissible};
use crate::expr::Expr;
use std::rc::Rc;

/// τ_Reg transformation: convert LCFL terms to left-linear form over tensor semiring
/// 
/// Given a term: a ⊗ Y ⊗ b (where Y is the differentiated variable)
/// Transform to: Z ⊗T (aᵗ ⊗ b) — left-linear in Z
/// 
/// The key insight: (aᵗ ⊗ b) is a tensor constant that "remembers" both coefficients.
/// Later, detensor-transpose recovers a ⊗ Y ⊗ b via (t,·)(aᵗ ⊗ b) when Y is substituted.

/// A term in the differentiated expression, representing c_left ⊗ Y ⊗ c_right
#[derive(Clone, Debug)]
pub struct LcflTerm<S: Semiring> {
    pub left: Expr<S>,   // coefficient on left of Y
    pub right: Expr<S>,  // coefficient on right of Y
}

impl<S: Semiring> LcflTerm<S> {
    pub fn new(left: Expr<S>, right: Expr<S>) -> Self {
        LcflTerm { left, right }
    }
    
    /// Identity term: 1 ⊗ Y ⊗ 1
    pub fn identity() -> Self {
        LcflTerm {
            left: Expr::one(),
            right: Expr::one(),
        }
    }
}

/// Computes the differential D_Xj[expr] and returns the LCFL coefficients directly.
/// 
/// This avoids building intermediate ASTs and completely prevents variable shadowing.
/// The key insight: when differentiating, the original variables become constants
/// (evaluated at ν), and only the linear placeholder Y appears in the result.
/// 
/// By computing coefficients directly, we:
/// - Eliminate the Star panic (D[a*] = a* ⊗ D[a] ⊗ a* wraps coefficients cleanly)
/// - Eliminate X vs Y shadowing (no AST to parse)
/// - Run much faster (no intermediate allocation)
pub fn differentiate_to_lcfl<S: Semiring>(
    expr: &Expr<S>,
    var_idx: usize
) -> Vec<LcflTerm<S>> {
    match expr {
        // Constants vanish in differentiation
        Expr::Zero | Expr::One | Expr::Const(_) => vec![],

        Expr::Var(i) => {
            if *i == var_idx {
                // D_X[X](Y) = Y = 1 ⊗ Y ⊗ 1
                vec![LcflTerm::identity()]
            } else {
                // Other variables act as constants, derivative is 0
                vec![]
            }
        }

        Expr::Combine(a, b) => {
            // D[a ⊕ b] = D[a] ⊕ D[b]
            let mut terms = differentiate_to_lcfl(a, var_idx);
            terms.extend(differentiate_to_lcfl(b, var_idx));
            terms
        }

        Expr::Extend(a, b) => {
            // Product rule: D[a ⊗ b] = D[a] ⊗ b ⊕ a ⊗ D[b]
            let mut terms = Vec::new();

            // For D[a] ⊗ b: the right coefficient gets multiplied by b
            for t in differentiate_to_lcfl(a, var_idx) {
                terms.push(LcflTerm::new(
                    t.left,
                    t.right.extend((**b).clone())
                ));
            }

            // For a ⊗ D[b]: the left coefficient gets multiplied by a
            for t in differentiate_to_lcfl(b, var_idx) {
                terms.push(LcflTerm::new(
                    (**a).clone().extend(t.left),
                    t.right
                ));
            }

            terms
        }

        Expr::Star(a) => {
            // Theorem 6.3: D[a*] = a* ⊗ D[a] ⊗ a*
            // This elegantly wraps the coefficients without panicking!
            let a_star = Expr::Star(a.clone());
            let mut terms = Vec::new();

            for t in differentiate_to_lcfl(a, var_idx) {
                terms.push(LcflTerm::new(
                    a_star.clone().extend(t.left),
                    t.right.extend(a_star.clone())
                ));
            }

            terms
        }
    }
}

/// Extract LCFL terms from a differentiated expression
/// 
/// After differentiation, we have expressions where Y appears linearly.
/// This function extracts all the (left, right) coefficient pairs.
/// 
/// Returns None if Y doesn't appear (constant term), Some(terms) otherwise.
pub fn extract_lcfl_terms<S: Semiring>(
    expr: &Expr<S>,
    var_idx: usize,
) -> Option<Vec<LcflTerm<S>>> {
    match expr {
        Expr::Zero => None, // No Y here
        Expr::One => None,  // No Y here (identity constant)
        Expr::Const(_) => None, // No Y here
        
        Expr::Var(i) => {
            if *i == var_idx {
                // Just Y by itself: 1 ⊗ Y ⊗ 1
                Some(vec![LcflTerm::identity()])
            } else {
                None // Different variable, treat as constant
            }
        }
        
        Expr::Combine(a, b) => {
            // Y appears in sum: collect from both sides
            let a_terms = extract_lcfl_terms(a, var_idx);
            let b_terms = extract_lcfl_terms(b, var_idx);
            
            match (a_terms, b_terms) {
                (None, None) => None,
                (Some(t), None) => Some(t),
                (None, Some(t)) => Some(t),
                (Some(mut a_t), Some(b_t)) => {
                    a_t.extend(b_t);
                    Some(a_t)
                }
            }
        }
        
        Expr::Extend(a, b) => {
            // Product: Y could be in a (giving c⊗Y on left) or b (giving Y⊗c on right)
            let a_terms = extract_lcfl_terms(a, var_idx);
            let b_terms = extract_lcfl_terms(b, var_idx);
            
            match (a_terms, b_terms) {
                (None, None) => None, // Y in neither - constant
                
                (Some(terms), None) => {
                    // Y in left: (c_l ⊗ Y ⊗ c_r) ⊗ b = c_l ⊗ Y ⊗ (c_r ⊗ b)
                    Some(terms.into_iter().map(|t| {
                        LcflTerm::new(
                            t.left,
                            t.right.extend((**b).clone()),
                        )
                    }).collect())
                }
                
                (None, Some(terms)) => {
                    // Y in right: a ⊗ (c_l ⊗ Y ⊗ c_r) = (a ⊗ c_l) ⊗ Y ⊗ c_r
                    Some(terms.into_iter().map(|t| {
                        LcflTerm::new(
                            (**a).clone().extend(t.left),
                            t.right,
                        )
                    }).collect())
                }
                
                (Some(_), Some(_)) => {
                    // Y appears in both sides - shouldn't happen for linear expressions
                    // This would be quadratic in Y
                    panic!("Non-linear expression: Y appears on both sides of extend")
                }
            }
        }
        
        Expr::Star(inner) => {
            // Star of something containing Y
            // From differentiation: g* ⊗ D[g] ⊗ g*
            // The D[g] part contains Y linearly
            // But the whole expression is still linear in Y
            
            // For now, we handle the case where the star itself doesn't contain Y
            // (the star gets evaluated before we see it in the differential)
            if extract_lcfl_terms(inner, var_idx).is_some() {
                // Y inside a star - this is the D[g*] = g* ⊗ D[g] ⊗ g* case
                // The structure is already linear, but nested
                // We'd need to recursively extract and compose
                panic!("Y inside Star - need to handle D[g*] case specially")
            } else {
                None // Star of constant
            }
        }
    }
}

/// Convert LCFL terms to tensor coefficients for left-linear system
/// 
/// Each term c_l ⊗ Y ⊗ c_r becomes coefficient (c_lᵗ ⊗ c_r) for Z
pub fn lcfl_to_tensor_coeffs<S: Admissible>(
    terms: &[LcflTerm<S>],
    values: &[S],
    one: &S,
) -> Vec<S::Tensor> {
    terms.iter().map(|term| {
        let left_val = term.left.eval_with_one(values, one);
        let right_val = term.right.eval_with_one(values, one);
        // Coupling: c_lᵗ ⊗ c_r
        left_val.transpose().tensor(&right_val)
    }).collect()
}

/// Sum tensor coefficients
pub fn sum_tensor_coeffs<S: Admissible>(coeffs: Vec<S::Tensor>) -> S::Tensor {
    coeffs.into_iter().fold(
        S::Tensor::zero(),
        |acc, c| acc.combine(&c)
    )
}

/// Full τ_Reg transformation for a differentiated equation
/// 
/// Input: D_Xj[Rhs_k]|ν(Y) — differentiated RHS for procedure k w.r.t. variable j
/// Output: Tensor coefficient for Z_j in the equation for Z_k
/// 
/// The full left-linear equation is:
///   Z_k = (1ᵗ ⊗ Rhs_k(ν)) ⊕T Σ_j (Z_j ⊗T Coeff_j(k))
/// 
/// where Coeff_j(k) = sum of (c_lᵗ ⊗ c_r) for all terms c_l ⊗ Y_j ⊗ c_r in D_Xj[Rhs_k]
/// Compute tau_reg coefficient directly from an expression (not pre-differentiated)
/// This uses the new differentiate_to_lcfl which avoids the X vs Y shadowing bug
/// `one` is a properly-sized identity element for the semiring
pub fn tau_reg_coefficient_direct<S: Admissible>(
    expr: &Expr<S>,
    var_idx: usize,
    values: &[S],
    one: &S,
) -> S::Tensor {
    let terms = differentiate_to_lcfl(expr, var_idx);
    if terms.is_empty() {
        S::Tensor::zero()
    } else {
        let coeffs = lcfl_to_tensor_coeffs(&terms, values, one);
        sum_tensor_coeffs::<S>(coeffs)
    }
}

/// Old interface for compatibility (takes pre-differentiated expression)
pub fn tau_reg_coefficient<S: Admissible>(
    diff_expr: &Expr<S>,
    var_idx: usize,
    values: &[S],
    one: &S,
) -> S::Tensor {
    match extract_lcfl_terms(diff_expr, var_idx) {
        None => S::Tensor::zero(),
        Some(terms) => {
            let coeffs = lcfl_to_tensor_coeffs(&terms, values, one);
            sum_tensor_coeffs::<S>(coeffs)
        }
    }
}

/// Compute the constant term (1ᵗ ⊗ Rhs(ν)) for the left-linear equation
/// Note: needs a "sized" identity, so we use rhs_value as the reference
pub fn tau_reg_constant<S: Admissible>(rhs_value: &S, one: &S) -> S::Tensor {
    one.transpose().tensor(rhs_value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean_matrix::BoolMatrix;
    use crate::differentiate::differentiate;
    
    #[test]
    fn test_extract_simple() {
        // Expression: X0 (just the variable)
        let x0 = Expr::<BoolMatrix>::var(0);
        let terms = extract_lcfl_terms(&x0, 0);
        
        assert!(terms.is_some());
        let terms = terms.unwrap();
        assert_eq!(terms.len(), 1);
        // Should be identity term: 1 ⊗ Y ⊗ 1
    }
    
    #[test]
    fn test_extract_product() {
        // Expression: c ⊗ X0 where c is a constant
        let c = Expr::<BoolMatrix>::constant(BoolMatrix::identity(2));
        let x0 = Expr::var(0);
        let expr = c.extend(x0);
        
        let terms = extract_lcfl_terms(&expr, 0);
        assert!(terms.is_some());
        let terms = terms.unwrap();
        assert_eq!(terms.len(), 1);
        // Should be: c ⊗ Y ⊗ 1
    }
    
    #[test]
    fn test_extract_sandwich() {
        // Expression: c ⊗ X0 ⊗ d (sandwich)
        let c = Expr::<BoolMatrix>::constant(BoolMatrix::identity(2));
        let d = Expr::<BoolMatrix>::constant(BoolMatrix::identity(2));
        let x0 = Expr::var(0);
        let expr = c.extend(x0).extend(d);
        
        let terms = extract_lcfl_terms(&expr, 0);
        assert!(terms.is_some());
        let terms = terms.unwrap();
        assert_eq!(terms.len(), 1);
        // Should be: c ⊗ Y ⊗ d
    }
    
    #[test]
    fn test_extract_sum() {
        // Expression: X0 ⊕ X0 (sum of same variable)
        let x0 = Expr::<BoolMatrix>::var(0);
        let expr = x0.clone().combine(x0);
        
        let terms = extract_lcfl_terms(&expr, 0);
        assert!(terms.is_some());
        let terms = terms.unwrap();
        assert_eq!(terms.len(), 2);
    }
    
    #[test]
    fn test_from_differentiation() {
        // Differentiate X0 ⊗ X0, then extract
        // D[X0 ⊗ X0] = 1 ⊗ X0 ⊕ X0 ⊗ 1
        let x0 = Expr::<BoolMatrix>::var(0);
        let expr = x0.clone().extend(x0);
        let diff = differentiate(&expr, 0);
        
        let terms = extract_lcfl_terms(&diff, 0);
        assert!(terms.is_some());
        let terms = terms.unwrap();
        // Should have 2 terms from product rule
        assert_eq!(terms.len(), 2);
    }
}
