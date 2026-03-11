use crate::semiring::{Semiring, Admissible};
use crate::boolean_matrix::BoolMatrix;
use crate::sparse_matrix::SparseBoolMatrix;
use crate::lifted_matrix::LiftedMatrix;
use crate::expr::Expr;
use crate::npa::solve_npa;

/// Simple LCG-based RNG for fuzzing (no external deps)
pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        // LCG parameters from Numerical Recipes
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.state
    }

    pub fn next_usize(&mut self, max: usize) -> usize {
        (self.next_u64() as usize) % max
    }

    pub fn next_bool(&mut self) -> bool {
        self.next_u64() & 1 == 1
    }

    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() as f64) / (u64::MAX as f64)
    }
}

/// Generate a random boolean matrix
pub fn random_bool_matrix(rng: &mut Rng, n: usize, density: f64) -> BoolMatrix {
    let mut m = BoolMatrix::new(n);
    for i in 0..n {
        for j in 0..n {
            if rng.next_f64() < density {
                m.set(i, j, true);
            }
        }
    }
    m
}

/// Generate a random sparse boolean matrix
pub fn random_sparse_matrix(rng: &mut Rng, n: usize, density: f64) -> SparseBoolMatrix {
    let mut m = SparseBoolMatrix::new(n);
    for i in 0..n {
        for j in 0..n {
            if rng.next_f64() < density {
                m.set(i, j, true);
            }
        }
    }
    m
}

/// Generate a random lifted matrix (with occasional Zero/One)
pub fn random_lifted_matrix(rng: &mut Rng, n: usize, density: f64) -> LiftedMatrix {
    match rng.next_usize(10) {
        0 => LiftedMatrix::Zero,
        1 => LiftedMatrix::One,
        _ => LiftedMatrix::Sized(random_sparse_matrix(rng, n, density)),
    }
}

/// Generate a random lifted matrix that's always Sized (for testing with abstract elements)
pub fn random_sized_lifted_matrix(rng: &mut Rng, n: usize, density: f64) -> LiftedMatrix {
    LiftedMatrix::Sized(random_sparse_matrix(rng, n, density))
}

/// Generate a random expression over n variables with given depth
pub fn random_expr<S: Semiring + Clone>(
    rng: &mut Rng,
    n_vars: usize,
    constants: &[S],
    max_depth: usize,
) -> Expr<S> {
    random_expr_impl(rng, n_vars, constants, max_depth, 0)
}

fn random_expr_impl<S: Semiring + Clone>(
    rng: &mut Rng,
    n_vars: usize,
    constants: &[S],
    max_depth: usize,
    depth: usize,
) -> Expr<S> {
    if depth >= max_depth {
        // At max depth, only generate leaves
        if n_vars > 0 && rng.next_bool() {
            Expr::var(rng.next_usize(n_vars))
        } else if !constants.is_empty() {
            Expr::constant(constants[rng.next_usize(constants.len())].clone())
        } else {
            Expr::one()
        }
    } else {
        match rng.next_usize(5) {
            0 => {
                // Variable
                if n_vars > 0 {
                    Expr::var(rng.next_usize(n_vars))
                } else {
                    Expr::one()
                }
            }
            1 => {
                // Constant
                if !constants.is_empty() {
                    Expr::constant(constants[rng.next_usize(constants.len())].clone())
                } else {
                    Expr::one()
                }
            }
            2 => {
                // Combine
                let a = random_expr_impl(rng, n_vars, constants, max_depth, depth + 1);
                let b = random_expr_impl(rng, n_vars, constants, max_depth, depth + 1);
                a.combine(b)
            }
            3 => {
                // Extend
                let a = random_expr_impl(rng, n_vars, constants, max_depth, depth + 1);
                let b = random_expr_impl(rng, n_vars, constants, max_depth, depth + 1);
                a.extend(b)
            }
            4 => {
                // Star (less common to avoid explosion)
                if rng.next_bool() {
                    let a = random_expr_impl(rng, n_vars, constants, max_depth, depth + 1);
                    a.star()
                } else {
                    random_expr_impl(rng, n_vars, constants, max_depth, depth + 1)
                }
            }
            _ => unreachable!()
        }
    }
}

/// Test semiring laws for a given semiring
pub fn test_semiring_laws<S: Semiring + std::fmt::Debug>(
    a: &S,
    b: &S, 
    c: &S,
) -> Result<(), String> {
    // Associativity of combine: (a ⊕ b) ⊕ c = a ⊕ (b ⊕ c)
    let lhs = a.combine(b).combine(c);
    let rhs = a.combine(&b.combine(c));
    if lhs != rhs {
        return Err(format!("combine not associative: {:?}", (a, b, c)));
    }

    // Commutativity of combine: a ⊕ b = b ⊕ a
    let lhs = a.combine(b);
    let rhs = b.combine(a);
    if lhs != rhs {
        return Err(format!("combine not commutative: {:?}", (a, b)));
    }

    // Identity for combine: a ⊕ 0 = a
    let lhs = a.combine(&S::zero());
    if lhs != *a {
        return Err(format!("zero not identity for combine: {:?}", a));
    }

    // Associativity of extend: (a ⊗ b) ⊗ c = a ⊗ (b ⊗ c)
    let lhs = a.extend(b).extend(c);
    let rhs = a.extend(&b.extend(c));
    if lhs != rhs {
        return Err(format!("extend not associative: {:?}", (a, b, c)));
    }

    // Identity for extend: a ⊗ 1 = a, 1 ⊗ a = a
    let one = S::one();
    let lhs = a.extend(&one);
    if lhs != *a {
        return Err(format!("one not right identity for extend: {:?}", a));
    }
    let lhs = one.extend(a);
    if lhs != *a {
        return Err(format!("one not left identity for extend: {:?}", a));
    }

    // Zero annihilates: a ⊗ 0 = 0, 0 ⊗ a = 0
    let zero = S::zero();
    let lhs = a.extend(&zero);
    if lhs != zero {
        return Err(format!("zero not right annihilator: {:?}", a));
    }
    let lhs = zero.extend(a);
    if lhs != zero {
        return Err(format!("zero not left annihilator: {:?}", a));
    }

    // Distributivity: a ⊗ (b ⊕ c) = (a ⊗ b) ⊕ (a ⊗ c)
    let lhs = a.extend(&b.combine(c));
    let rhs = a.extend(b).combine(&a.extend(c));
    if lhs != rhs {
        return Err(format!("extend not left-distributive: {:?}", (a, b, c)));
    }

    // Right distributivity: (a ⊕ b) ⊗ c = (a ⊗ c) ⊕ (b ⊗ c)
    let lhs = a.combine(b).extend(c);
    let rhs = a.extend(c).combine(&b.extend(c));
    if lhs != rhs {
        return Err(format!("extend not right-distributive: {:?}", (a, b, c)));
    }

    Ok(())
}

/// Test admissible semiring laws
pub fn test_admissible_laws<S: Admissible + std::fmt::Debug>(
    a: &S,
    b: &S,
) -> Result<(), String> {
    // Transpose is involution: (aᵗ)ᵗ = a
    let lhs = a.transpose().transpose();
    if lhs != *a {
        return Err(format!("transpose not involution: {:?}", a));
    }

    // Transpose distributes over combine: (a ⊕ b)ᵗ = aᵗ ⊕ bᵗ
    let lhs = a.combine(b).transpose();
    let rhs = a.transpose().combine(&b.transpose());
    if lhs != rhs {
        return Err(format!("transpose doesn't distribute over combine: {:?}", (a, b)));
    }

    // Transpose reverses extend: (a ⊗ b)ᵗ = bᵗ ⊗ aᵗ
    let lhs = a.extend(b).transpose();
    let rhs = b.transpose().extend(&a.transpose());
    if lhs != rhs {
        return Err(format!("transpose doesn't reverse extend: {:?}", (a, b)));
    }

    // Detensor-transpose recovers extend: (t,·)(aᵗ ⊗ b) = a ⊗ b
    // This is the fundamental theorem of NPA-TP (Equation 43)
    let tensor_prod = a.transpose().tensor(b);
    let recovered = S::detensor_transpose(&tensor_prod);
    let expected = a.extend(b);
    if recovered != expected {
        return Err(format!(
            "detensor-transpose failed: (t,·)(aᵗ ⊗ b) = {:?} != a ⊗ b = {:?}",
            recovered, expected
        ));
    }

    Ok(())
}

/// Test that sparse and dense matrices produce equivalent results
pub fn test_sparse_dense_equivalence(
    rng: &mut Rng,
    n: usize,
    density: f64,
) -> Result<(), String> {
    let dense_a = random_bool_matrix(rng, n, density);
    let dense_b = random_bool_matrix(rng, n, density);

    // Convert to sparse
    let sparse_a = dense_to_sparse(&dense_a);
    let sparse_b = dense_to_sparse(&dense_b);

    // Test combine
    let dense_c = dense_a.combine(&dense_b);
    let sparse_c = sparse_a.combine(&sparse_b);
    if dense_to_sparse(&dense_c) != sparse_c {
        return Err("combine mismatch".to_string());
    }

    // Test extend
    let dense_c = dense_a.extend(&dense_b);
    let sparse_c = sparse_a.extend(&sparse_b);
    if dense_to_sparse(&dense_c) != sparse_c {
        return Err("extend mismatch".to_string());
    }

    // Test transpose
    let dense_c = dense_a.transpose();
    let sparse_c = sparse_a.transpose();
    if dense_to_sparse(&dense_c) != sparse_c {
        return Err("transpose mismatch".to_string());
    }

    // Test star
    let dense_c = dense_a.star();
    let sparse_c = sparse_a.star();
    if dense_to_sparse(&dense_c) != sparse_c {
        return Err("star mismatch".to_string());
    }

    Ok(())
}

fn dense_to_sparse(d: &BoolMatrix) -> SparseBoolMatrix {
    let mut s = SparseBoolMatrix::new(d.n);
    for i in 0..d.n {
        for j in 0..d.n {
            if d.get(i, j) {
                s.set(i, j, true);
            }
        }
    }
    s
}

/// Generate a LINEAR expression suitable for NPA
/// Each variable appears at most once, and not inside stars
pub fn random_linear_expr<S: Semiring + Clone>(
    rng: &mut Rng,
    n_vars: usize,
    constants: &[S],
    max_depth: usize,
) -> Expr<S> {
    let mut used_vars = std::collections::HashSet::default();
    random_linear_expr_impl(rng, n_vars, constants, max_depth, 0, &mut used_vars)
}

fn random_linear_expr_impl<S: Semiring + Clone>(
    rng: &mut Rng,
    n_vars: usize,
    constants: &[S],
    max_depth: usize,
    depth: usize,
    used_vars: &mut std::collections::HashSet<usize>,
) -> Expr<S> {
    if depth >= max_depth {
        // At max depth, only constants
        if !constants.is_empty() {
            Expr::constant(constants[rng.next_usize(constants.len())].clone())
        } else {
            Expr::one()
        }
    } else {
        match rng.next_usize(5) {
            0 => {
                // Variable (only if not already used)
                let available: Vec<usize> = (0..n_vars)
                    .filter(|v| !used_vars.contains(v))
                    .collect();
                if !available.is_empty() && rng.next_bool() {
                    let var = available[rng.next_usize(available.len())];
                    used_vars.insert(var);
                    Expr::var(var)
                } else if !constants.is_empty() {
                    Expr::constant(constants[rng.next_usize(constants.len())].clone())
                } else {
                    Expr::one()
                }
            }
            1 => {
                // Constant
                if !constants.is_empty() {
                    Expr::constant(constants[rng.next_usize(constants.len())].clone())
                } else {
                    Expr::one()
                }
            }
            2 => {
                // Combine (linear in both branches)
                let a = random_linear_expr_impl(rng, n_vars, constants, max_depth, depth + 1, used_vars);
                let b = random_linear_expr_impl(rng, n_vars, constants, max_depth, depth + 1, used_vars);
                a.combine(b)
            }
            3 => {
                // Extend with variable on one side only (LCFL form: a ⊗ X ⊗ b)
                let available: Vec<usize> = (0..n_vars)
                    .filter(|v| !used_vars.contains(v))
                    .collect();
                
                if !available.is_empty() && rng.next_bool() {
                    let var = available[rng.next_usize(available.len())];
                    used_vars.insert(var);
                    
                    // Generate a ⊗ X ⊗ b (sandwich form)
                    let a = random_const_expr(rng, constants, max_depth, depth + 1);
                    let b = random_const_expr(rng, constants, max_depth, depth + 1);
                    a.extend(Expr::var(var)).extend(b)
                } else {
                    // Just constants
                    let a = random_const_expr(rng, constants, max_depth, depth + 1);
                    let b = random_const_expr(rng, constants, max_depth, depth + 1);
                    a.extend(b)
                }
            }
            4 => {
                // Star of CONSTANTS only (no variables inside star)
                let a = random_const_expr(rng, constants, max_depth, depth + 1);
                a.star()
            }
            _ => unreachable!()
        }
    }
}

/// Generate an expression with only constants (no variables)
fn random_const_expr<S: Semiring + Clone>(
    rng: &mut Rng,
    constants: &[S],
    max_depth: usize,
    depth: usize,
) -> Expr<S> {
    if depth >= max_depth || rng.next_usize(3) == 0 {
        if !constants.is_empty() {
            Expr::constant(constants[rng.next_usize(constants.len())].clone())
        } else {
            Expr::one()
        }
    } else {
        match rng.next_usize(3) {
            0 => {
                let a = random_const_expr(rng, constants, max_depth, depth + 1);
                let b = random_const_expr(rng, constants, max_depth, depth + 1);
                a.combine(b)
            }
            1 => {
                let a = random_const_expr(rng, constants, max_depth, depth + 1);
                let b = random_const_expr(rng, constants, max_depth, depth + 1);
                a.extend(b)
            }
            2 => {
                let a = random_const_expr(rng, constants, max_depth, depth + 1);
                a.star()
            }
            _ => unreachable!()
        }
    }
}

/// Naive Kleene iteration: repeatedly evaluate equations until fixpoint
/// This is trivially correct for finite lattices (Boolean matrices)
fn naive_kleene_iteration(
    rhs: &[Expr<BoolMatrix>],
    max_rounds: usize,
    mat_size: usize,
) -> Option<Vec<BoolMatrix>> {
    if rhs.is_empty() {
        return Some(vec![]);
    }
    
    let n = rhs.len();
    
    // Initialize to properly-sized zero matrices
    let mut values: Vec<BoolMatrix> = vec![BoolMatrix::new(mat_size); n];
    
    for _ in 0..max_rounds {
        let new_values: Vec<BoolMatrix> = rhs.iter()
            .map(|expr| expr.eval(&values))
            .collect();
        
        // Combine with previous (monotonic update)
        let combined: Vec<BoolMatrix> = values.iter()
            .zip(new_values.iter())
            .map(|(old, new)| old.combine(new))
            .collect();
        
        if combined == values {
            return Some(combined);
        }
        
        values = combined;
    }
    
    None // didn't converge
}

/// Test NPA convergence on random equation systems
pub fn test_npa_convergence(
    rng: &mut Rng,
    n_equations: usize,
    matrix_size: usize,
    expr_depth: usize,
    max_rounds: usize,
) -> Result<(), String> {
    // Generate random constants
    let constants: Vec<BoolMatrix> = (0..3)
        .map(|_| random_bool_matrix(rng, matrix_size, 0.3))
        .collect();

    // Generate random LINEAR RHS expressions (suitable for NPA)
    let rhs: Vec<Expr<BoolMatrix>> = (0..n_equations)
        .map(|_| random_linear_expr(rng, n_equations, &constants, expr_depth))
        .collect();

    let one = BoolMatrix::identity(matrix_size);
    
    // Run NPA
    let result = solve_npa(rhs, one, max_rounds);

    // Check it converged
    if result.rounds > max_rounds {
        return Err(format!("NPA didn't converge in {} rounds", max_rounds));
    }

    // Check values are valid (non-zero dimension)
    for (i, v) in result.values.iter().enumerate() {
        if v.n != matrix_size && v.n != 0 {
            return Err(format!("NPA result {} has wrong dimension: {}", i, v.n));
        }
    }

    Ok(())
}

/// Test semiring laws for LiftedMatrix (handles Zero/One properly)
pub fn test_lifted_semiring_laws(
    a: &LiftedMatrix,
    b: &LiftedMatrix, 
    c: &LiftedMatrix,
) -> Result<(), String> {
    // Associativity of combine: (a ⊕ b) ⊕ c = a ⊕ (b ⊕ c)
    let lhs = a.combine(b).combine(c);
    let rhs = a.combine(&b.combine(c));
    if lhs != rhs {
        return Err(format!("combine not associative"));
    }

    // Commutativity of combine: a ⊕ b = b ⊕ a
    let lhs = a.combine(b);
    let rhs = b.combine(a);
    if lhs != rhs {
        return Err(format!("combine not commutative"));
    }

    // Identity for combine: a ⊕ 0 = a
    let lhs = a.combine(&LiftedMatrix::Zero);
    if lhs != *a {
        return Err(format!("zero not identity for combine"));
    }

    // Associativity of extend: (a ⊗ b) ⊗ c = a ⊗ (b ⊗ c)
    let lhs = a.extend(b).extend(c);
    let rhs = a.extend(&b.extend(c));
    if lhs != rhs {
        return Err(format!("extend not associative"));
    }

    // Identity for extend: a ⊗ 1 = a, 1 ⊗ a = a
    let one = LiftedMatrix::One;
    let lhs = a.extend(&one);
    if lhs != *a {
        return Err(format!("one not right identity for extend"));
    }
    let lhs = one.extend(a);
    if lhs != *a {
        return Err(format!("one not left identity for extend"));
    }

    // Zero annihilates: a ⊗ 0 = 0, 0 ⊗ a = 0
    let zero = LiftedMatrix::Zero;
    let lhs = a.extend(&zero);
    if lhs != zero {
        return Err(format!("zero not right annihilator"));
    }
    let lhs = zero.extend(a);
    if lhs != zero {
        return Err(format!("zero not left annihilator"));
    }

    // Distributivity: a ⊗ (b ⊕ c) = (a ⊗ b) ⊕ (a ⊗ c)
    let lhs = a.extend(&b.combine(c));
    let rhs = a.extend(b).combine(&a.extend(c));
    if lhs != rhs {
        return Err(format!("extend not left-distributive"));
    }

    // Right distributivity: (a ⊕ b) ⊗ c = (a ⊗ c) ⊕ (b ⊗ c)
    let lhs = a.combine(b).extend(c);
    let rhs = a.extend(c).combine(&b.extend(c));
    if lhs != rhs {
        return Err(format!("extend not right-distributive"));
    }

    Ok(())
}

/// Run the full fuzzer suite
pub fn run_fuzz(seed: u64, iterations: usize) -> Vec<String> {
    let mut rng = Rng::new(seed);
    let mut failures = Vec::new();

    for i in 0..iterations {
        let n = 2 + rng.next_usize(3); // 2-4

        // Test LiftedMatrix semiring laws (with abstract Zero/One)
        let a = random_lifted_matrix(&mut rng, n, 0.3);
        let b = random_lifted_matrix(&mut rng, n, 0.3);
        let c = random_lifted_matrix(&mut rng, n, 0.3);

        if let Err(e) = test_lifted_semiring_laws(&a, &b, &c) {
            failures.push(format!("iter {}: lifted semiring law: {}", i, e));
        }

        // Test sparse/dense equivalence (using raw sparse matrices)
        if let Err(e) = test_sparse_dense_equivalence(&mut rng, n, 0.3) {
            failures.push(format!("iter {}: sparse/dense: {}", i, e));
        }

        // Test NPA convergence with LiftedMatrix (less frequently due to cost)
        // Note: NPA still uses BoolMatrix directly, so we test that separately
        if i % 10 == 0 {
            let n_eq = 1 + rng.next_usize(7); // test up to 8 equations
            let mat_size = 2;
            if let Err(e) = test_npa_convergence(&mut rng, n_eq, mat_size, 3, 50) {
                failures.push(format!("iter {}: NPA: {}", i, e));
            }
        }
    }

    failures
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rng() {
        let mut rng = Rng::new(42);
        let a = rng.next_u64();
        let b = rng.next_u64();
        assert_ne!(a, b);
    }

    #[test]
    fn test_random_matrix() {
        let mut rng = Rng::new(42);
        let m = random_bool_matrix(&mut rng, 3, 0.5);
        assert_eq!(m.n, 3);
    }

    #[test]
    fn test_random_expr() {
        let mut rng = Rng::new(42);
        let constants = vec![BoolMatrix::identity(2)];
        let expr: Expr<BoolMatrix> = random_expr(&mut rng, 2, &constants, 3);
        // Just check it doesn't panic
        let values = vec![BoolMatrix::identity(2), BoolMatrix::identity(2)];
        let _ = expr.eval(&values);
    }

    #[test]
    fn test_fuzz_quick() {
        let failures = run_fuzz(12345, 100);
        if !failures.is_empty() {
            for f in &failures {
                eprintln!("{}", f);
            }
            panic!("{} failures", failures.len());
        }
    }

    #[test]
    fn test_fuzz_extended() {
        let failures = run_fuzz(98765, 500);
        if !failures.is_empty() {
            for f in &failures {
                eprintln!("{}", f);
            }
            panic!("{} failures", failures.len());
        }
    }

    #[test]
    #[ignore] // run with: cargo test --release -- --ignored
    fn test_fuzz_long() {
        let failures = run_fuzz(31415926, 50000);
        if !failures.is_empty() {
            for f in &failures[..std::cmp::min(100, failures.len())] {
                eprintln!("{}", f);
            }
            panic!("{} failures", failures.len());
        }
    }

    #[test]
    #[ignore] // run with: cargo test --release -- --ignored test_npa_vs_naive
    fn test_npa_vs_naive() {
        use super::*;
        let mut rng = Rng::new(27182818);
        let mut failures = Vec::new();

        for i in 0..1000 {
            let n_eq = 1 + rng.next_usize(7); // test up to 8 equations
            let mat_size = 2 + rng.next_usize(14); // test up to 16x16 matrices
            
            let constants: Vec<BoolMatrix> = (0..3)
                .map(|_| random_bool_matrix(&mut rng, mat_size, 0.3))
                .collect();

            // Use random_expr (non-linear) to stress test NPA properly
            // NPA's power is handling quadratic equations and stars
            let rhs: Vec<Expr<BoolMatrix>> = (0..n_eq)
                .map(|_| random_expr(&mut rng, n_eq, &constants, 3))
                .collect();

            let one = BoolMatrix::identity(mat_size);
            let npa_result = solve_npa(rhs.clone(), one, 100);
            
            if let Some(naive_result) = naive_kleene_iteration(&rhs, 1000, mat_size) {
                if npa_result.values != naive_result {
                    failures.push(format!(
                        "iter {}: NPA != naive\n  NPA: {:?}\n  naive: {:?}",
                        i, npa_result.values, naive_result
                    ));
                }
            } else {
                failures.push(format!("iter {}: naive didn't converge", i));
            }
        }

        if !failures.is_empty() {
            for f in &failures[..std::cmp::min(10, failures.len())] {
                eprintln!("{}", f);
            }
            panic!("{} failures", failures.len());
        }
    }
}
