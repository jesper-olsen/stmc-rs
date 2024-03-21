use marsaglia_rs::marsaglia::Marsaglia;
use std::f64::consts::PI;
use marsaglia_rs::gamma::gamma_ln;

fn gamma_ln_test() {
    // ln(1), ln(1), ln(2) => 0.0, 0.0, 0.693 ...
    let x_values = vec![1.0, 2.0, 3.0];
    for x in x_values {
        println!("Gamma_ln({}) = {}", x, gamma_ln(x));
    }
}

// C MOMENTS FOR NORMALLY DISTRIBUTED RANDOM NUMBERS.
fn main() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    println!("\n CALCULATION OF NORMAL GAUSSIAN MOMENTS 1, 2, 4, 6 AND 8:");
    println!("\n K    NDAT      G1         G2        G4        G6        G8");
    for k in 1..=10 {
        let mut xmo = [0.0f64; 5];
        let mut xmom = [0.0f64; 5];

        let ndat = 2i32.pow(2 * k - 1);
        xmo[0] = 1.0;
        for _ in 0..ndat {
            let x = rng.gauss();
            xmom[0] += x;
            for i in 1..5 {
                xmo[i] = xmo[i - 1] * x.powi(2);
                xmom[i] += xmo[i];
            }
        }
        // normalisation
        for e in &mut xmom {
            *e /= ndat as f64
        }
        println!(
            "{:2} {:8} {:12.4} {:7.4} {:10.4} {:10.4} {:10.4}",
            k, ndat as i64, xmom[0], xmom[1], xmom[2], xmom[3], xmom[4]
        );
    }

    let mut xmom = [0.0f64; 5];
    for i in 1..5 {
        let xn = 2.0f64 * (i as f64);
        let g = gamma_ln(0.5 * xn + 0.5).exp();
        xmom[i] = 2.0f64.powf(0.5 * xn) * g / PI.sqrt();
    }
    println!(
        "     EXACT: {:12.4} {:7.4} {:10.4} {:10.4} {:10.4}",
        xmom[0], xmom[1], xmom[2], xmom[3], xmom[4]
    );
}
