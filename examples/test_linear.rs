// minimal test case
use dataflow::*;
use dataflow::linear::{extract_linear_system, solve_linear, is_linear, extract_linear};

fn main() {
    // Simple: X = c ⊕ a⊗X
    // Solution: X = a*⊗c
    
    let c = BoolMatrix::identity(2); // I
    let mut a = BoolMatrix::new(2);
    a.set(0, 1, true); // off-diagonal only
    
    println!("a = {:?}", a);
    println!("c = {:?}", c);
    println!("a* = {:?}", a.star());
    
    // Build: c ⊕ (a ⊗ X)
    let expr: Expr<BoolMatrix> = Expr::Combine(
        Expr::Const(c.clone()).into(),
        Expr::Extend(
            Expr::Const(a.clone()).into(),
            Expr::Var(0).into(),
        ).into(),
    );
    
    let rhs = vec![expr];
    
    // Extraction
    let (constants, coeffs) = extract_linear_system(&rhs);
    println!("\nExtracted:");
    println!("  constants[0] = {:?}", constants[0]);
    println!("  coeffs[0][0] = {:?}", coeffs[0][0]);
    
    // Expected: constants = [I], coeffs = [[a]]
    // A* should be [[a*]] = [[I ⊕ a]]
    // X = A*·c = (I ⊕ a)·I = I ⊕ a
    
    let one = BoolMatrix::identity(2);
    let linear = solve_linear(constants.clone(), coeffs.clone(), one.clone());
    println!("\nLinear solver result: {:?}", linear[0]);
    
    // NPA for comparison
    let npa = solve_npa(rhs.clone(), one, 100);
    println!("NPA result: {:?}", npa.values[0]);
    
    // Expected
    let expected = a.star().extend(&c);
    println!("Expected (a*·c): {:?}", expected);
}
