use dataflow::*;

fn main() {
    let mut a = BoolMatrix::new(2);
    a.set(0, 0, true);
    let b = BoolMatrix::identity(2);
    
    // X0 = a* ⊗ X0 ⊗ b
    let rhs0: Expr<BoolMatrix> = 
        Expr::Star(std::rc::Rc::new(Expr::constant(a.clone())))
            .extend(Expr::var(0))
            .extend(Expr::constant(b.clone()));
    
    let one = BoolMatrix::identity(2);
    
    let npa = solve_npa(vec![rhs0.clone()], one.clone(), 10);
    println!("NPA ({} rounds): {:?}", npa.rounds, npa.values[0]);
    
    // Naive
    let zero = BoolMatrix::new(2);
    let mut x = zero.clone();
    for round in 0..10 {
        let new_x = rhs0.eval(&[x.clone()]);
        let combined = x.combine(&new_x);
        if combined == x {
            println!("Naive converged at round {}", round);
            break;
        }
        x = combined;
    }
    println!("Naive: {:?}", x);
    
    println!("\nMatch: {}", npa.values[0] == x);
}
