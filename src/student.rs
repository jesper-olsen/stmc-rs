use crate::beta::{beta, beta_i};
use crate::chi2::fi1;

pub fn stud_df(t: f64, nf: usize) -> f64 {
    let a = 0.5 * nf as f64;
    let b = 0.5;
    let x = (nf as f64) / (nf as f64 + t.powi(2));
    match t {
        t if t < 0.0 => 0.5 * beta_i(x, a, b),
        t if t > 0.0 => 1.0 - 0.5 * beta_i(x, a, b),
        _ => 0.5,
    }
}

pub fn stud_pd(t: f64, nf: usize) -> f64 {
    let f = nf as f64;
    let fhalf = 0.5 * f;
    let sf = f.sqrt();
    (1.0 + (t / sf).powi(2)).powf(-(f + 1.0) / 2.0) / sf / beta(0.5, fhalf)
}

pub fn stud_qdf(t: f64, nf: usize) -> f64 {
    let f = stud_df(t, nf);
    if f > 0.5 {
        1.0 - f
    } else {
        f
    }
}

pub fn stud_xq(q: f64, nf: usize) -> f64 {
    if q == 0.5 {
        return 0.0;
    }

    let mut x1 = 0.0;
    let mut x2 = 0.0;

    if q > 0.5 {
        loop {
            x2 += 1.0;
            if stud_df(x2, nf) > q {
                break;
            }
        }
    } else {
        loop {
            x1 -= 1.0;
            if stud_df(x1, nf) < q {
                break;
            }
        }
    }

    let f = |x: f64| stud_df(x, nf);
    fi1(f, q, x1, x2)
}
