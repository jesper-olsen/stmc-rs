const ITER_MAX: usize = 800;
const EPS: f64 = 3.0e-9;

// C LN OF GAMMA FUNCTION ALA LANCZOS, SIAM Num. Anal. B1 (1964) 1.
pub fn gamma_ln(x: f64) -> f64 {
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

pub fn gamma_p(a: f64, x: f64) -> f64 {
    if x < 0.0 || a <= 0.0 {
        println!("A,X: {}, {}", a, x);
        panic!("GAMMA_P - INPUT NOT COVERED");
    }

    let gln = gamma_ln(a);
    let mut lbad1 = false;
    let mut lbad2 = false;
    loop {
        if x < a + 1.0 && !lbad2 || lbad1 {
            // series expansion
            if x < 0.0 {
                return 0.0;
            }

            let mut sum = 1.0 / a;
            let mut add = sum;
            for i in 1..=ITER_MAX {
                add *= x / (a + i as f64);
                sum += add;
                if add.abs() < sum.abs() * EPS {
                    return sum * (-x + a * x.ln() - gln).exp();
                }
            }
            if lbad1 {
                panic!("GAMMA_P - 1: A too large, ITER_MAX too small")
            }
            lbad2 = true;
        } else {
            // continued fraction expansion
            let mut gold = 0.0;
            let mut a0 = x + 1.0 - a;
            let mut b0 = 1.0;
            let mut a1 = x * a0 + x;
            let mut b1 = x * b0 + 1.0;
            for i in 2..=ITER_MAX {
                let xiter = i as f64;
                let xiterma = xiter - a;
                a0 = (a1 + a0 * xiterma) / a1;
                b0 = (b1 + b0 * xiterma) / a1;
                b1 = x * b0 + xiter * b1 / a1;
                a1 = x * a0 + xiter;
                if a1 == 0.0 {
                    panic!("a1 zero in gamma_p");
                }
                let g = b1 / a1;
                if ((g - gold) / g).abs() < EPS {
                    return 1.0 - (-x + a * x.ln() - gln).exp() * g; // and take its complement
                }
                gold = g;
            }

            if lbad2 {
                panic!("GAMMA_P - 2: A too large, ITER_MAX too small")
            }

            lbad1 = true;
        }
    }
}

pub fn error_f(x: f64) -> f64 {
    let g = gamma_p(0.5, x.powi(2));
    match x {
        x if x > 0.0 => g,
        x if x < 0.0 => -g,
        _ => 0.0,
    }
}
