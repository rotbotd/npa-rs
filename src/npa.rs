use crate::semiring::{Semiring, Admissible};
use crate::expr::Expr;
use crate::differentiate::differentiate;
use crate::regularize::{tau_reg_coefficient, tau_reg_constant};

/// Result of NPA-TP analysis
#[derive(Clone, Debug)]
pub struct NpaResult<S: Semiring> {
    pub values: Vec<S>,
    pub rounds: usize,
}

/// NPA-TP solver for interprocedural dataflow analysis
/// 
/// Implements Algorithm 7.1 from the POPL 2016 paper.
pub struct NpaSolver<S: Admissible> {
    /// Number of procedures/equations
    pub n: usize,
    /// Right-hand side expressions for each equation
    pub rhs: Vec<Expr<S>>,
    /// Path expressions from Tarjan (precomputed)
    pub path_exprs: Vec<Expr<S::Tensor>>,
    /// Identity element (properly sized)
    pub one: S,
}

impl<S: Admissible> NpaSolver<S> {
    /// Create a new solver from a system of equations
    /// 
    /// Each equation has form: X_i = rhs_i(X_0, ..., X_{n-1})
    /// `one` must be a properly-sized identity element
    pub fn new(rhs: Vec<Expr<S>>, one: S) -> Self {
        let n = rhs.len();
        
        // Build the dependence graph for the linearized system
        // This is a simplified version - full version would do Tarjan on the dependence graph
        let path_exprs = Self::build_path_expressions(n, &rhs);
        
        NpaSolver { n, rhs, path_exprs, one }
    }

    /// Build parameterized path expressions for the tensor system
    /// 
    /// In the full algorithm, we'd:
    /// 1. Build dependence graph G with edges Z_k → Z_j labeled <k,j>
    /// 2. Run Tarjan on G to get regular expressions R_i
    /// 
    /// For now, we use a simplified direct iteration approach.
    fn build_path_expressions(n: usize, _rhs: &[Expr<S>]) -> Vec<Expr<S::Tensor>> {
        // Simplified: just return placeholders
        // The real implementation would construct the dependence graph
        // and run Tarjan to get the parameterized expressions
        vec![Expr::zero(); n]
    }

    /// Run Newton iteration until convergence
    pub fn solve(&self, max_rounds: usize) -> NpaResult<S> {
        // Initialize: ν = f(0)
        let zero_values: Vec<S> = vec![S::zero(); self.n];
        let mut nu: Vec<S> = self.rhs.iter()
            .map(|rhs| rhs.eval(&zero_values))
            .collect();
        
        let mut round = 0;
        
        loop {
            round += 1;
            if round > max_rounds {
                break;
            }
            
            // Compute tensor coefficients T_kj at current ν
            let coeffs = self.compute_coefficients(&nu);
            
            // Solve left-linear system in tensor semiring
            let z_values = self.solve_tensor_system(&nu, &coeffs);
            
            // Apply detensor-transpose to get new ν
            let new_nu: Vec<S> = z_values.iter()
                .map(|z| S::detensor_transpose(z))
                .collect();
            
            // Check convergence
            if new_nu == nu {
                nu = new_nu;
                break;
            }
            
            nu = new_nu;
        }
        
        NpaResult { values: nu, rounds: round }
    }

    /// Compute tensor coefficients T_kj = Coeff_k(τReg(D_Xk[Rhs_j]|ν))
    fn compute_coefficients(&self, nu: &[S]) -> Vec<Vec<S::Tensor>> {
        let mut coeffs = vec![vec![S::Tensor::zero(); self.n]; self.n];
        
        for j in 0..self.n {
            for k in 0..self.n {
                // Differentiate Rhs_j with respect to X_k
                let diff = differentiate(&self.rhs[j], k);
                // Apply τ_Reg to get tensor coefficient
                coeffs[k][j] = tau_reg_coefficient(&diff, k, nu);
            }
        }
        
        coeffs
    }

    /// Solve the left-linear system over tensor semiring
    /// 
    /// Z_j = (1ᵗ ⊗ Rhs_j(ν)) ⊕T Σ_k (Z_k ⊗T T_kj)
    /// 
    /// For now, use simple iteration. The full algorithm uses Tarjan's
    /// precomputed path expressions.
    fn solve_tensor_system(
        &self, 
        nu: &[S], 
        coeffs: &[Vec<S::Tensor>],
    ) -> Vec<S::Tensor> {
        // Evaluate RHS at current nu
        let rhs_values: Vec<S> = self.rhs.iter().map(|rhs| rhs.eval(nu)).collect();
        
        // Initialize with constant terms using properly sized identity
        let mut z: Vec<S::Tensor> = rhs_values.iter()
            .map(|rhs_val| tau_reg_constant(rhs_val, &self.one))
            .collect();
        
        // Iterate until fixpoint (simplified Kleene iteration in tensor space)
        let mut prev = vec![S::Tensor::zero(); self.n];
        
        while z != prev {
            prev = z.clone();
            
            for j in 0..self.n {
                // Z_j = constant ⊕T Σ_k (Z_k ⊗T T_kj)
                let mut new_zj = tau_reg_constant(&rhs_values[j], &self.one);
                
                for k in 0..self.n {
                    let term = prev[k].extend(&coeffs[k][j]);
                    new_zj = new_zj.combine(&term);
                }
                
                z[j] = new_zj;
            }
        }
        
        z
    }
}

/// Convenience function to solve a system of equations
/// `one` must be a properly-sized identity element
pub fn solve_npa<S: Admissible>(rhs: Vec<Expr<S>>, one: S, max_rounds: usize) -> NpaResult<S> {
    let solver = NpaSolver::new(rhs, one);
    solver.solve(max_rounds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean_matrix::BoolMatrix;
    
    #[test]
    fn test_simple_equation() {
        // X = a (constant equation, no recursion)
        let a = BoolMatrix::identity(2);
        let one = BoolMatrix::identity(2);
        let rhs = vec![Expr::constant(a.clone())];
        
        let result = solve_npa(rhs, one, 10);
        
        assert_eq!(result.values.len(), 1);
        assert_eq!(result.values[0], a);
        // Constant equation converges in 1 round
    }
    
    #[test]
    fn test_self_reference() {
        // X = a ⊕ X (should converge to a*)
        // But since X appears linearly, Newton should handle it
        let a = BoolMatrix::identity(2);
        let one = BoolMatrix::identity(2);
        let rhs = vec![
            Expr::constant(a.clone()).combine(Expr::var(0))
        ];
        
        let result = solve_npa(rhs, one, 20);
        
        // Should converge
        assert!(result.rounds <= 20);
    }
    
    #[test]
    fn test_two_equations() {
        // X0 = a
        // X1 = X0 ⊗ b
        let a = BoolMatrix::identity(2);
        let b = BoolMatrix::identity(2);
        let one = BoolMatrix::identity(2);
        
        let rhs = vec![
            Expr::constant(a.clone()),
            Expr::var(0).extend(Expr::constant(b.clone())),
        ];
        
        let result = solve_npa(rhs, one, 10);
        
        assert_eq!(result.values.len(), 2);
        // X0 = a = identity
        assert_eq!(result.values[0], a);
        // X1 = X0 ⊗ b = identity ⊗ identity = identity
        assert_eq!(result.values[1], a.extend(&b));
    }
}

#[cfg(test)]
mod debug_tests {
    use super::*;
    use crate::boolean_matrix::BoolMatrix;
    use crate::regularize::tau_reg_constant;
    
    #[test]
    fn debug_simple() {
        let a = BoolMatrix::identity(2);
        let one = BoolMatrix::identity(2);
        println!("a = {:?}", a);
        
        // tau_reg_constant creates 1ᵗ ⊗ a
        let tensor = one.transpose().tensor(&a);
        println!("tensor (manual) = {:?}", tensor);
        
        // tau_reg_constant now takes sized one
        let tensor2 = tau_reg_constant(&a, &one);
        println!("tensor (tau_reg) = {:?}", tensor2);
        
        // detensor should recover a
        let recovered = BoolMatrix::detensor_transpose(&tensor);
        println!("recovered = {:?}", recovered);
        
        assert_eq!(recovered, a);
        assert_eq!(tensor, tensor2);
    }
}
