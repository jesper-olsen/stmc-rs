use marsaglia_rs::Marsaglia;
use std::f64::consts::PI;

// C LN OF GAMMA FUNCTION ALA LANCZOS, SIAM Num. Anal. B1 (1964) 1.
fn gamma_ln(x: f64) -> f64 {
    const C1_L: f64 = 76.18009173;
    const C2_L: f64 = -86.50532033;
    const C3_L: f64 = 24.01409822;
    const C4_L: f64 = -1.231739516;
    const C5_L: f64 = 0.120858003e-2;
    const C6_L: f64 = -0.536382e-5;
    const STP_L: f64 = 2.50662827465;

    if x <= 0.0 {
        panic!("GAMMA_LN: Argument X = {}", x);
    }
    let y = if x > 1.0 {
        x // Full accuracy of Lanczos formula.
    } else {
        1.0 + x // Use Gamma(z+1)=z*Gamma(z).
    };
    let ser = ((1.0 + C1_L / y) + C2_L / (y + 1.0)) + C3_L / (y + 2.0);
    let ser = ((ser + C4_L / (y + 3.0)) + C5_L / (y + 4.0)) + C6_L / (y + 5.0);
    let gamma_ln = (y - 0.5) * (y + 4.5).ln() - (y + 4.5) + (STP_L * ser).ln();
    if x > 1.0 {
        gamma_ln
    } else {
        gamma_ln - x.ln()
    }
}

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
