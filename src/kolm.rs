use crate::gamma::gamma_ln;

// Exact one-sided Kolmogorov tests, implementing equation of
// Birnbaum and Tingey, Ann. Math. Stat. 22 (1951) 592. See also:
// van der Waerden, Mathematical Statistics, Springer 1969.

pub fn kolm1(fxct: &[f64]) -> (f64, f64, f64, f64) {
    let mut del1 = 0.0;
    let mut del2 = 0.0;
    let mut q1 = 0.0;
    let mut q2 = 0.0;
    let n = fxct.len();
    for j in 1..=2 {
        let mut del = 0.0f64;
        if j == 1 {
            for i in 0..n {
                let femp = i as f64 / n as f64;
                del = del.max(fxct[i] - femp);
            }
            del1 = del;
        } else {
            for i in 0..n {
                let femp = (i + 1) as f64 / n as f64;
                del = del.max(femp - fxct[i]);
            }
            del2 = del;
        }

        let delln = del.ln();
        let xn = (n + 1) as f64;
        let ynln = gamma_ln(xn);
        let mut q = 0.0;
        for k in 0..=n {
            let xkp1 = (k + 1) as f64;
            let xnp1mk = (n + 1 - k) as f64;
            let bln = ynln - gamma_ln(xkp1) - gamma_ln(xnp1mk);
            let delplus = del + k as f64 / n as f64;
            let x = 1.0 - delplus;
            if x >= 0.0 {
                let qln =
                    bln + delln + (k as i32 - 1) as f64 * delplus.ln() + (n - k) as f64 * x.ln();
                q += qln.exp();
            }
        }

        if j == 1 {
            q1 = q
        } else {
            q2 = q;
        }
    }
    (del1, del2, q1, q2)
}

// Asymptotic two-sided Kolmogorov test in the form of
// M.A. Stephens, J. Royal Stat. Soc. B 32 (1970) 115.
pub fn kolm2_as(fxct: &[f64]) -> (f64, f64) {
    let mut del = 0.0f64;
    let n = fxct.len();
    for i in 0..n {
        let femp = i as f64 / n as f64;
        del = del.max(fxct[i] - femp);
        let femp = (i + 1) as f64 / n as f64;
        del = del.max(femp - fxct[i]);
    }
    let sqn = (n as f64).sqrt();
    let a = -2.0 * (sqn * del + (12.0 * del) / 100.0 + (11.0 * del) / (100.0 * sqn)).powi(2);
    let mut sign_two = 2.0;
    let mut q = 0.0f64;
    let mut cut = 0.0;
    for j in 1..=100 {
        let add = sign_two * (a * (j as i32).pow(2) as f64).exp();
        q += add;
        if add.abs() < cut {
            return (del, q);
        }
        sign_two = -sign_two;
        cut = add.abs() / 1000.0;
    }
    (del, 1.0)
}
