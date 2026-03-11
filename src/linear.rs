//! Linear equation solver using Gaussian elimination
//!
//! For systems where each variable appears at most once per RHS (no X*X terms),
//! we can solve directly without Newton iteration.
//!
//! System form: X_i = c_i ⊕ (a_i1 ⊗ X_1) ⊕ (a_i2 ⊗ X_2) ⊕ ... ⊕ (a_in ⊗ X_n)
//!
//! In matrix form: X = c ⊕ A·X
//! Solution: X = A* · c  (where A* is Kleene star / transitive closure)

use crate::expr::Expr;
use crate::semiring::Semiring;

/// Check if an expression is linear in the given variables
/// (each variable appears at most once, no stars containing variables)
pub fn is_linear<S: Semiring>(expr: &Expr<S>, n_vars: usize) -> bool {
    let mut var_counts = vec![0usize; n_vars];
    is_linear_inner(expr, &mut var_counts)
}

/// Check if an expression is LEFT-linear: variables only appear as `coeff ⊗ X`,
/// never as `X ⊗ coeff` or `c1 ⊗ X ⊗ c2`.
/// 
/// This is required for extract_linear() and solve_linear().
pub fn is_left_linear<S: Semiring>(expr: &Expr<S>, n_vars: usize) -> bool {
    is_left_linear_inner(expr, n_vars, true)
}

fn is_left_linear_inner<S: Semiring>(expr: &Expr<S>, n_vars: usize, can_have_var: bool) -> bool {
    match expr {
        Expr::Zero | Expr::One | Expr::Const(_) => true,
        Expr::Var(i) => can_have_var && *i < n_vars,
        Expr::Combine(a, b) => {
            is_left_linear_inner(a, n_vars, can_have_var) && 
            is_left_linear_inner(b, n_vars, can_have_var)
        }
        Expr::Extend(a, b) => {
            // For left-linear: a must be variable-free, b can have variables
            // This ensures variables only appear on the RIGHT of extends
            let a_has_var = has_variable(a, n_vars);
            if a_has_var {
                // Variable on left side of extend - not left-linear
                false
            } else {
                // a is constant, recurse into b
                is_left_linear_inner(b, n_vars, can_have_var)
            }
        }
        Expr::Star(inner) => {
            // Star must be variable-free for linearity
            !has_variable(inner, n_vars)
        }
    }
}

fn has_variable<S: Semiring>(expr: &Expr<S>, n_vars: usize) -> bool {
    match expr {
        Expr::Zero | Expr::One | Expr::Const(_) => false,
        Expr::Var(i) => *i < n_vars,
        Expr::Combine(a, b) | Expr::Extend(a, b) => {
            has_variable(a, n_vars) || has_variable(b, n_vars)
        }
        Expr::Star(inner) => has_variable(inner, n_vars),
    }
}

fn is_linear_inner<S: Semiring>(expr: &Expr<S>, var_counts: &mut [usize]) -> bool {
    match expr {
        Expr::Zero | Expr::One | Expr::Const(_) => true,
        Expr::Var(i) => {
            if *i < var_counts.len() {
                var_counts[*i] += 1;
                var_counts[*i] <= 1
            } else {
                false // out of bounds variable
            }
        }
        Expr::Combine(a, b) => {
            is_linear_inner(a, var_counts) && is_linear_inner(b, var_counts)
        }
        Expr::Extend(a, b) => {
            // Product is linear only if at most one side has variables
            let mut a_counts = var_counts.to_vec();
            let mut b_counts = var_counts.to_vec();
            
            if !is_linear_inner(a, &mut a_counts) || !is_linear_inner(b, &mut b_counts) {
                return false;
            }
            
            // Check that variables from a and b don't overlap (would make X*X)
            for i in 0..var_counts.len() {
                let a_has = a_counts[i] > var_counts[i];
                let b_has = b_counts[i] > var_counts[i];
                if a_has && b_has {
                    return false; // same var on both sides = quadratic
                }
                var_counts[i] = a_counts[i].max(b_counts[i]);
            }
            true
        }
        Expr::Star(inner) => {
            // Star is only linear if it contains no variables
            let mut inner_counts = var_counts.to_vec();
            if !is_linear_inner(inner, &mut inner_counts) {
                return false;
            }
            // Check no new variables were found
            for i in 0..var_counts.len() {
                if inner_counts[i] > var_counts[i] {
                    return false; // variable inside star = nonlinear
                }
            }
            true
        }
    }
}

/// Result of extracting a left-linear expression
#[derive(Clone, Debug)]
pub struct LinearForm<S: Semiring> {
    pub constant: S,
    pub coeffs: Vec<S>,
}

/// Extract linear form from a LEFT-linear expression.
/// 
/// For an expression of the form `c ⊕ Σ(a_j ⊗ X_j)`, returns:
/// - constant = c
/// - coeffs[j] = a_j
///
/// IMPORTANT: Only works for left-linear expressions where variables appear
/// as `coeff ⊗ X`, never `X ⊗ coeff`. Use is_left_linear() to check first.
/// Sandwich forms like `c1 ⊗ X ⊗ c2` will give incorrect results.
///
/// For non-left-linear systems, use NPA instead.
pub fn extract_linear<S: Semiring + Clone>(expr: &Expr<S>, n_vars: usize) -> LinearForm<S> {
    extract_linear_with_coeff(expr, n_vars, S::one())
}

fn extract_linear_with_coeff<S: Semiring + Clone>(
    expr: &Expr<S>,
    n_vars: usize,
    left_coeff: S,
) -> LinearForm<S> {
    match expr {
        Expr::Zero => LinearForm {
            constant: S::zero(),
            coeffs: vec![S::zero(); n_vars],
        },
        Expr::One => LinearForm {
            constant: left_coeff,
            coeffs: vec![S::zero(); n_vars],
        },
        Expr::Const(c) => LinearForm {
            constant: left_coeff.extend(c),
            coeffs: vec![S::zero(); n_vars],
        },
        Expr::Var(j) => {
            let mut coeffs = vec![S::zero(); n_vars];
            if *j < n_vars {
                coeffs[*j] = left_coeff;
            }
            LinearForm {
                constant: S::zero(),
                coeffs,
            }
        }
        Expr::Combine(a, b) => {
            let la = extract_linear_with_coeff(a, n_vars, left_coeff.clone());
            let lb = extract_linear_with_coeff(b, n_vars, left_coeff);
            LinearForm {
                constant: la.constant.combine(&lb.constant),
                coeffs: la.coeffs.iter().zip(lb.coeffs.iter())
                    .map(|(a, b)| a.combine(b))
                    .collect(),
            }
        }
        Expr::Extend(a, b) => {
            // For left-linear forms: only b can contain variables.
            // a ⊗ b with a constant, b = (c ⊕ Σ a_j⊗X_j)
            // = a⊗c ⊕ Σ (a⊗a_j)⊗X_j
            //
            // So we recurse on b with left_coeff multiplied by a's value.
            let a_const = eval_constant_part(a);
            extract_linear_with_coeff(b, n_vars, left_coeff.extend(&a_const))
        }
        Expr::Star(inner) => {
            // For linear expressions, star can only contain constants
            let inner_const = eval_constant_part(inner);
            LinearForm {
                constant: left_coeff.extend(&inner_const.star()),
                coeffs: vec![S::zero(); n_vars],
            }
        }
    }
}

/// Extract the constant term and coefficient matrix from linear expressions
/// 
/// For X_i = c_i ⊕ Σ(a_ij ⊗ X_j), returns (c, A) where:
/// - c[i] = constant term for equation i  
/// - A[i][j] = coefficient of X_j in equation i
pub fn extract_linear_system<S: Semiring + Clone>(
    rhs: &[Expr<S>],
) -> (Vec<S>, Vec<Vec<S>>) {
    let n = rhs.len();
    let mut constants = vec![S::zero(); n];
    let mut coeffs = vec![vec![S::zero(); n]; n];
    
    for i in 0..n {
        let form = extract_linear(&rhs[i], n);
        constants[i] = form.constant;
        for j in 0..n {
            coeffs[i][j] = form.coeffs[j].clone();
        }
    }
    
    (constants, coeffs)
}

/// Evaluate an expression treating all variables as zero
fn eval_constant_part<S: Semiring + Clone>(expr: &Expr<S>) -> S {
    match expr {
        Expr::Zero => S::zero(),
        Expr::One => S::one(),
        Expr::Const(c) => c.clone(),
        Expr::Var(_) => S::zero(),
        Expr::Combine(a, b) => eval_constant_part(a).combine(&eval_constant_part(b)),
        Expr::Extend(a, b) => eval_constant_part(a).extend(&eval_constant_part(b)),
        Expr::Star(inner) => eval_constant_part(inner).star(),
    }
}

/// Solve a linear system X = c ⊕ A·X
/// 
/// Solution: X = A* · c
/// 
/// We compute A* using iteration, then multiply by c.
pub fn solve_linear<S: Semiring + Clone + PartialEq>(
    constants: Vec<S>,
    coeffs: Vec<Vec<S>>,
    one: S, // identity element of proper size
) -> Vec<S> {
    let n = constants.len();
    if n == 0 {
        return vec![];
    }
    
    // Compute A* (transitive closure of coefficient matrix)
    let a_star = matrix_star(coeffs, one);
    
    // X = A* · c
    // result[i] = Σ_j (A*[i][j] ⊗ c[j])
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        // Start with first term, then combine rest
        let mut acc = a_star[i][0].extend(&constants[0]);
        for j in 1..n {
            acc = acc.combine(&a_star[i][j].extend(&constants[j]));
        }
        result.push(acc);
    }
    
    result
}

/// Compute matrix Kleene star (transitive closure) using iteration
/// A* = I ⊕ A ⊕ A² ⊕ A³ ⊕ ...
/// 
/// `one` is the identity element of proper size (e.g., 2×2 identity matrix)
fn matrix_star<S: Semiring + Clone + PartialEq>(a: Vec<Vec<S>>, one: S) -> Vec<Vec<S>> {
    let n = a.len();
    if n == 0 {
        return vec![];
    }
    
    // Initialize result = I ⊕ A
    // Diagonal entries get `one ⊕ a[i][i]`, off-diagonal get `a[i][j]`
    let mut result: Vec<Vec<S>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(n);
        for j in 0..n {
            if i == j {
                // Diagonal: I ⊕ A[i][i] = one ⊕ a[i][i]
                row.push(one.clone().combine(&a[i][j]));
            } else {
                // Off-diagonal: 0 ⊕ A[i][j] = a[i][j]
                row.push(a[i][j].clone());
            }
        }
        result.push(row);
    }
    
    // Iterate until fixpoint
    for _ in 0..n*n*10 {
        let prev = result.clone();
        
        // result = result ⊕ result·A
        let ra = matrix_mult(&result, &a);
        for i in 0..n {
            for j in 0..n {
                result[i][j] = result[i][j].combine(&ra[i][j]);
            }
        }
        
        // Check convergence
        if result == prev {
            break;
        }
    }
    
    result
}

fn matrix_mult<S: Semiring + Clone>(a: &[Vec<S>], b: &[Vec<S>]) -> Vec<Vec<S>> {
    let n = a.len();
    let mut result = vec![vec![S::zero(); n]; n];
    
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                result[i][j] = result[i][j].combine(&a[i][k].extend(&b[k][j]));
            }
        }
    }
    
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BoolMatrix;

    #[test]
    fn test_is_linear_simple() {
        // X = a (constant) - linear
        let a = BoolMatrix::identity(2);
        let expr: Expr<BoolMatrix> = Expr::Const(a);
        assert!(is_linear(&expr, 1));
    }

    #[test]
    fn test_is_linear_with_var() {
        // X = a ⊗ Y - linear
        let a = BoolMatrix::identity(2);
        let expr: Expr<BoolMatrix> = Expr::Extend(
            Expr::Const(a).into(),
            Expr::Var(0).into(),
        );
        assert!(is_linear(&expr, 1));
    }

    #[test]
    fn test_is_linear_quadratic() {
        // X = Y ⊗ Y - NOT linear
        let expr: Expr<BoolMatrix> = Expr::Extend(
            Expr::Var(0).into(),
            Expr::Var(0).into(),
        );
        assert!(!is_linear(&expr, 1));
    }

    #[test]
    fn test_is_linear_star_of_var() {
        // X = Y* - NOT linear
        let expr: Expr<BoolMatrix> = Expr::Star(Expr::Var(0).into());
        assert!(!is_linear(&expr, 1));
    }

    #[test]
    fn test_is_linear_star_of_const() {
        // X = a* - linear (star of constant is fine)
        let a = BoolMatrix::identity(2);
        let expr: Expr<BoolMatrix> = Expr::Star(Expr::Const(a).into());
        assert!(is_linear(&expr, 1));
    }

    #[test]
    fn test_sandwich_not_left_linear() {
        // c1 ⊗ X ⊗ c2 - variable in the middle, NOT left-linear
        let mut c1 = BoolMatrix::new(2);
        c1.set(0, 0, true);
        let mut c2 = BoolMatrix::new(2);
        c2.set(1, 1, true);
        
        // Build: c1 ⊗ (X ⊗ c2) = Extend(c1, Extend(X, c2))
        let sandwich: Expr<BoolMatrix> = Expr::Extend(
            Expr::Const(c1.clone()).into(),
            Expr::Extend(
                Expr::Var(0).into(),
                Expr::Const(c2.clone()).into(),
            ).into(),
        );
        
        // This IS linear (each var appears once, no vars in stars)
        assert!(is_linear(&sandwich, 1));
        
        // But it's NOT left-linear (var has stuff on its right)
        assert!(!is_left_linear(&sandwich, 1));
        
        // Left-linear form: c1 ⊗ X
        let left_linear: Expr<BoolMatrix> = Expr::Extend(
            Expr::Const(c1.clone()).into(),
            Expr::Var(0).into(),
        );
        assert!(is_left_linear(&left_linear, 1));
    }

    #[test]
    fn test_extract_linear_simple() {
        // X = a ⊕ b⊗Y
        // constant = a, coeff[0] = b
        let a = BoolMatrix::identity(2);
        let mut b = BoolMatrix::new(2);
        b.set(0, 0, true);
        b.set(1, 1, true);
        
        let expr: Expr<BoolMatrix> = Expr::Combine(
            Expr::Const(a.clone()).into(),
            Expr::Extend(
                Expr::Const(b.clone()).into(),
                Expr::Var(0).into(),
            ).into(),
        );
        
        let form = extract_linear(&expr, 1);
        assert_eq!(form.constant, a);
        assert_eq!(form.coeffs[0], b);
    }

    #[test]
    fn test_solve_linear_vs_npa() {
        use crate::npa::solve_npa;
        use crate::fuzz::{Rng, random_bool_matrix};
        
        let mut rng = Rng::new(12345);
        
        for _ in 0..100 {
            let mat_size = 2 + rng.next_usize(4);
            let n_eq = 1 + rng.next_usize(4);
            
            // Generate random LINEAR system: X_i = c_i ⊕ Σ(a_ij ⊗ X_j)
            let constants: Vec<BoolMatrix> = (0..n_eq)
                .map(|_| random_bool_matrix(&mut rng, mat_size, 0.3))
                .collect();
            let coeffs: Vec<Vec<BoolMatrix>> = (0..n_eq)
                .map(|_| (0..n_eq)
                    .map(|_| random_bool_matrix(&mut rng, mat_size, 0.2))
                    .collect())
                .collect();
            
            // Build expressions
            let rhs: Vec<Expr<BoolMatrix>> = (0..n_eq)
                .map(|i| {
                    let mut expr = Expr::Const(constants[i].clone());
                    for j in 0..n_eq {
                        expr = Expr::Combine(
                            expr.into(),
                            Expr::Extend(
                                Expr::Const(coeffs[i][j].clone()).into(),
                                Expr::Var(j).into(),
                            ).into(),
                        );
                    }
                    expr
                })
                .collect();
            
            // Verify all are linear
            for expr in &rhs {
                assert!(is_linear(expr, n_eq), "Generated non-linear expression");
            }
            
            // Solve with linear solver
            let one = BoolMatrix::identity(mat_size);
            let (lin_const, lin_coeffs) = extract_linear_system(&rhs);
            let linear_result = solve_linear(lin_const, lin_coeffs, one.clone());
            
            // Solve with NPA
            let npa_result = solve_npa(rhs, one, 100);
            
            // Compare
            assert_eq!(linear_result, npa_result.values, 
                "Linear solver disagrees with NPA");
        }
    }
}
