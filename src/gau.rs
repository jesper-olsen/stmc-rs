use crate::chi2::chi2pdf_xq;
use crate::chi2::fi1;
use crate::gamma;
use std::f64::consts::PI;

/// COMPARISION OF TWO MEANS: (GAUSSIAN DIFFERENCE TEST).
/// INPUT: TWO GAUSSIAN DATA POINTS, MEAN VALUES AND ERROR BARS.
/// OUPUT: LIKELIHOOD Q  THAT THE DISCREPANCY IS DUE TO CHANCE.
pub fn gaudif(xm1: f64, eb1: f64, xm2: f64, eb2: f64) -> f64 {
    let sigma = (eb1.powi(2) + eb2.powi(2)).sqrt();
    let xx = (xm1 - xm2).abs() / (sigma * 2.0f64.sqrt());
    1.0 - gamma::error_f(xx)
}

/// GAUSSIAN PROBABILITY DENSITY FUNCTION.
pub fn gau_pd(x: f64, mean: f64, std_dev: f64) -> f64 {
    (-0.5 * (x - mean).powi(2)).exp() / (std_dev * (2.0 * PI).sqrt())
}

/// GAUSSIAN CUMULATIVE DISTRIBUTION FUNCTION.
pub fn gau_df(x: f64) -> f64 {
    //0.5 *(1.0 + erf(x / (2.0f64.sqrt())))
    0.5 * (1.0 + gamma::error_f(x / (2.0f64.sqrt())))
}

/// GAUSSIAN PEAKED DISTRIBUTION FUNCTION.
pub fn gau_qdf(x: f64) -> f64 {
    let f = 0.5 + 0.5 * gamma::error_f(x / 2.0f64.sqrt());
    if f <= 0.5 {
        f
    } else {
        1.0 - f
    }
}

pub fn gau_xq(q: f64) -> f64 {
    if q == 0.5 {
        return 0.0;
    }

    let mut x1 = 0.0;
    let mut x2 = 0.0;

    if q > 0.5 {
        loop {
            x2 += 1.0;
            if gau_df(x2) > q {
                break;
            }
        }
    } else {
        loop {
            x1 -= 1.0;
            if gau_df(x1) < q {
                break;
            }
        }
    }

    fi1(gau_df, q, x1, x2)
}

/// ERROR BAR FOR THE GAUSSIAN ERROR BAR.
/// INPUT:  NUMBER OF GAUSSIAN DATA, CONFIDENCE INTERVAL.
/// OUTPUT: UPPER AND LOWER ERROR BAR LIMITS (ASYMPTOTIC FORMULA).
pub fn sebar_e_as(n: usize, pc: f64) -> (f64, f64) {
    let sdv = if pc <= 0.0 {
        gau_xq(0.5 * 1.68268949)
    } else {
        gau_xq(0.5 * (1.0 + pc))
    };
    let d = sdv / (2.0 * (n - 1) as f64).sqrt();
    (1.0 + d, 1.0 - d)
}

//ERROR BAR FOR THE GAUSSIAN ERROR BAR.
//INPUT:   NUMBER OF GAUSSIAN DATA, CONFIDENCE INTERVAL.
//OUTPUT:  UPPER AND LOWER ERROR BAR LIMITS (EXACT).
pub fn sebar_e(n: usize, pc: f64) -> (f64, f64) {
    let nf = n - 1;
    if n > 17000 {
        sebar_e_as(n, pc)
    } else {
        let q = if pc <= 0.0 {
            0.5 * (1.0 - 0.954499736)
        } else {
            0.5 * (1.0 - pc)
        };
        let p = 1.0 - q;
        (
            1.0 / chi2pdf_xq(q, nf).sqrt(),
            1.0 / chi2pdf_xq(p, nf).sqrt(),
        )
    }
}
