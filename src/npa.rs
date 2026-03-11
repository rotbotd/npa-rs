use rustc_hash::FxHashMap as HashMap;
use crate::semiring::{Semiring, Admissible};
use crate::expr::Expr;
use crate::cfg::{Cfg, DomTree};
use crate::tarjan::tarjan;
use crate::differentiate::differentiate;
use crate::regularize::{differentiate_to_lcfl, lcfl_to_tensor_coeffs, sum_tensor_coeffs, tau_reg_constant, LcflTerm};

/// Result of NPA-TP analysis
#[derive(Clone, Debug)]
pub struct NpaResult<S: Semiring> {
    pub values: Vec<S>,
    pub rounds: usize,
}

/// Label for dependence graph edges: <k, j> means Z_k appears in equation for Z_j
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DepLabel {
    pub k: usize,  // source variable
    pub j: usize,  // target equation
}

/// DepLabel as a "fake" semiring — we never actually evaluate these,
/// just build expressions over them and substitute later
impl Semiring for DepLabel {
    fn zero() -> Self {
        DepLabel { k: usize::MAX, j: usize::MAX }
    }
    
    fn one() -> Self {
        DepLabel { k: usize::MAX - 1, j: usize::MAX - 1 }
    }
    
    fn combine(&self, _other: &Self) -> Self {
        panic!("DepLabel::combine should never be called")
    }
    
    fn extend(&self, _other: &Self) -> Self {
        panic!("DepLabel::extend should never be called")  
    }
    
    fn star(&self) -> Self {
        panic!("DepLabel::star should never be called")
    }
}

/// NPA-TP solver for interprocedural dataflow analysis
/// 
/// Implements Algorithm 7.1 from the POPL 2016 paper.
pub struct NpaSolver<S: Admissible> {
    /// Number of procedures/equations
    pub n: usize,
    /// Right-hand side expressions for each equation
    pub rhs: Vec<Expr<S>>,
    /// Dependence graph (computed once)
    pub dep_graph: Cfg<DepLabel>,
    /// Path expressions from Tarjan on dependence graph (computed once)
    pub path_exprs: HashMap<usize, Expr<DepLabel>>,
    /// Precomputed LCFL terms for each (k, j) coefficient
    /// lcfl_terms[k][j] = differentiate_to_lcfl(rhs[j], k)
    pub lcfl_terms: Vec<Vec<Vec<LcflTerm<S>>>>,
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
        
        // Build dependence graph G:
        // - Dummy vertex 0 (entry)
        // - Vertices 1..=n for Z_1..Z_n
        // - Edge from 0 to j labeled <0,j> for constant term
        // - Edge from k to j labeled <k,j> if Z_k appears in RHS of Z_j
        let dep_graph = Self::build_dependence_graph(n, &rhs);
        
        // Run Tarjan on dependence graph to get path expressions
        let domtree = DomTree::compute(&dep_graph);
        let path_exprs = tarjan(&dep_graph, &domtree);
        
        // Precompute LCFL terms for all (k, j) coefficient pairs
        // The derivative structure never changes, only the values we evaluate at
        let mut lcfl_terms = vec![vec![vec![]; n]; n];
        for j in 0..n {
            for k in 0..n {
                lcfl_terms[k][j] = differentiate_to_lcfl(&rhs[j], k);
            }
        }
        
        NpaSolver { n, rhs, dep_graph, path_exprs, lcfl_terms, one }
    }

    /// Build dependence graph from RHS expressions
    fn build_dependence_graph(n: usize, rhs: &[Expr<S>]) -> Cfg<DepLabel> {
        // Vertex 0 = dummy entry (Λ in the paper)
        // Vertices 1..=n = Z_1..Z_n
        let mut cfg = Cfg::new(0);
        
        for j in 1..=n {
            cfg.add_node(j);
            
            // Edge from entry to each Z_j (constant term)
            cfg.add_edge(0, j, DepLabel { k: 0, j });
            
            // Check which variables appear in RHS of equation j-1
            let vars = collect_vars(&rhs[j - 1]);
            for k in vars {
                // k is 0-indexed variable, add 1 for vertex numbering
                cfg.add_edge(k + 1, j, DepLabel { k: k + 1, j });
            }
        }
        
        cfg
    }

    /// Run Newton iteration until convergence
    pub fn solve(&self, max_rounds: usize) -> NpaResult<S> {
        // Initialize: ν = f(0)
        // Use properly-sized zeros, not sentinel zeros
        let size = self.one.size();
        let zero_values: Vec<S> = vec![S::zero_sized(size); self.n];
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
            
            // Solve using path expressions from Tarjan
            let z_values = self.solve_with_path_exprs(&nu, &coeffs);
            
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
                // Use precomputed LCFL terms - only evaluate at current nu
                let terms = &self.lcfl_terms[k][j];
                if terms.is_empty() {
                    coeffs[k][j] = S::Tensor::zero();
                } else {
                    let tensor_coeffs = lcfl_to_tensor_coeffs(terms, nu, &self.one);
                    coeffs[k][j] = sum_tensor_coeffs::<S>(tensor_coeffs);
                }
            }
        }
        
        coeffs
    }

    /// Solve using precomputed path expressions from Tarjan
    /// 
    /// Each Z_j = R_j[<0,j> ← (1ᵗ ⊗ Rhs_j(ν)), <k,j> ← T_kj]
    fn solve_with_path_exprs(
        &self, 
        nu: &[S], 
        coeffs: &[Vec<S::Tensor>],
    ) -> Vec<S::Tensor> {
        let rhs_values: Vec<S> = self.rhs.iter().map(|rhs| rhs.eval(nu)).collect();
        
        // Build substitution map for path expression evaluation
        let mut label_values: HashMap<DepLabel, S::Tensor> = HashMap::default();
        
        for j in 1..=self.n {
            // Constant term: <0,j> → (1ᵗ ⊗ Rhs_j(ν))
            label_values.insert(
                DepLabel { k: 0, j },
                tau_reg_constant(&rhs_values[j - 1], &self.one),
            );
            
            // Variable terms: <k,j> → T_{k-1,j-1}
            for k in 1..=self.n {
                label_values.insert(
                    DepLabel { k, j },
                    coeffs[k - 1][j - 1].clone(),
                );
            }
        }
        
        // Evaluate path expressions for each Z_j
        let mut z_values = Vec::with_capacity(self.n);
        
        // Get tensor size from the one matrix (it knows its dimension)
        let tensor_size = self.one.tensor(&self.one).size();
        
        for j in 1..=self.n {
            if let Some(path_expr) = self.path_exprs.get(&j) {
                let z_j = eval_path_expr::<S>(path_expr, &label_values, tensor_size);
                z_values.push(z_j);
            } else {
                // No path to this node - use constant term only
                z_values.push(tau_reg_constant(&rhs_values[j - 1], &self.one));
            }
        }
        
        z_values
    }

    /// Debug: print path expressions
    pub fn debug_path_exprs(&self) {
        eprintln!("Path expressions:");
        for j in 1..=self.n {
            if let Some(path_expr) = self.path_exprs.get(&j) {
                eprintln!("  Z{} = {:?}", j, path_expr);
            } else {
                eprintln!("  Z{} = (no path)", j);
            }
        }
    }
}

/// Collect all variable indices that appear in an expression
fn collect_vars<S: Semiring>(expr: &Expr<S>) -> Vec<usize> {
    let mut vars = Vec::new();
    collect_vars_impl(expr, &mut vars);
    vars.sort();
    vars.dedup();
    vars
}

fn collect_vars_impl<S: Semiring>(expr: &Expr<S>, vars: &mut Vec<usize>) {
    match expr {
        Expr::Zero | Expr::One | Expr::Const(_) => {}
        Expr::Var(i) => vars.push(*i),
        Expr::Combine(a, b) | Expr::Extend(a, b) => {
            collect_vars_impl(a, vars);
            collect_vars_impl(b, vars);
        }
        Expr::Star(a) => collect_vars_impl(a, vars),
    }
}

/// Evaluate a path expression over DepLabel alphabet with given substitution
/// `tensor_size` is the size of the tensor matrices (n^2 where n is the base matrix size)
fn eval_path_expr<S: Admissible>(
    expr: &Expr<DepLabel>,
    label_values: &HashMap<DepLabel, S::Tensor>,
    tensor_size: usize,
) -> S::Tensor {
    // Use memoization to avoid O(2^N) blowup on DAG expressions
    let mut memo: HashMap<*const Expr<DepLabel>, S::Tensor> = HashMap::default();
    eval_path_expr_memo::<S>(expr, label_values, tensor_size, &mut memo)
}

fn eval_path_expr_memo<S: Admissible>(
    expr: &Expr<DepLabel>,
    label_values: &HashMap<DepLabel, S::Tensor>,
    tensor_size: usize,
    memo: &mut HashMap<*const Expr<DepLabel>, S::Tensor>,
) -> S::Tensor {
    // Check memo first using pointer address as key
    let ptr = expr as *const Expr<DepLabel>;
    if let Some(cached) = memo.get(&ptr) {
        return cached.clone();
    }
    
    let result = match expr {
        Expr::Zero => S::Tensor::zero_sized(tensor_size),
        Expr::One => S::Tensor::one_sized(tensor_size),
        Expr::Const(label) => {
            // Check for special sentinel values from DepLabel semiring
            if *label == DepLabel::zero() {
                S::Tensor::zero_sized(tensor_size)
            } else if *label == DepLabel::one() {
                S::Tensor::one_sized(tensor_size)
            } else {
                label_values.get(label).cloned().unwrap_or_else(|| S::Tensor::zero_sized(tensor_size))
            }
        }
        Expr::Var(_) => {
            // Variables in path expressions shouldn't exist after Tarjan
            S::Tensor::zero_sized(tensor_size)
        }
        Expr::Combine(a, b) => {
            let a_val: S::Tensor = eval_path_expr_memo::<S>(a, label_values, tensor_size, memo);
            let b_val: S::Tensor = eval_path_expr_memo::<S>(b, label_values, tensor_size, memo);
            a_val.combine(&b_val)
        }
        Expr::Extend(a, b) => {
            let a_val: S::Tensor = eval_path_expr_memo::<S>(a, label_values, tensor_size, memo);
            let b_val: S::Tensor = eval_path_expr_memo::<S>(b, label_values, tensor_size, memo);
            a_val.extend(&b_val)
        }
        Expr::Star(a) => {
            let a_val: S::Tensor = eval_path_expr_memo::<S>(a, label_values, tensor_size, memo);
            a_val.star()
        }
    };
    
    memo.insert(ptr, result.clone());
    result
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
    }
    
    #[test]
    fn test_self_reference() {
        // X = a ⊕ X (should converge to a*)
        let a = BoolMatrix::identity(2);
        let one = BoolMatrix::identity(2);
        let rhs = vec![
            Expr::constant(a.clone()).combine(Expr::var(0))
        ];
        
        let result = solve_npa(rhs, one, 20);
        
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
        assert_eq!(result.values[0], a);
        assert_eq!(result.values[1], a.extend(&b));
    }
    
    #[test]
    fn test_collect_vars() {
        let expr: Expr<BoolMatrix> = Expr::var(0)
            .extend(Expr::var(1))
            .combine(Expr::var(0));
        
        let vars = collect_vars(&expr);
        assert_eq!(vars, vec![0, 1]);
    }
}
