use crate::gamma::{gamma_ln, gamma_p};

pub fn chi2_df(chi2: f64, nf: usize) -> f64 {
    let a = 0.5 * nf as f64;
    let x = 0.5 * chi2;
    gamma_p(a, x)
}

pub fn chi2pdf_df(chi2: f64, nf: usize) -> f64 {
    let a = 0.5 * nf as f64;
    let x = a * chi2;
    gamma_p(a, x)
}

pub fn chi2pdf_pd(chi2: f64, nf: usize) -> f64 {
    let a = 0.5 * nf as f64;
    let aln = a.ln();
    let cln = chi2.ln();
    let yln = aln - a * chi2 + (a - 1.0) * (aln + cln) - gamma_ln(a);
    yln.exp()
}

/// peaked distribution function.
pub fn chi2pdf_qdf(chi2: f64, nf: usize) -> f64 {
    let a = 0.5 * nf as f64;
    let x = a * chi2;
    let g = gamma_p(a, x);
    if g <= 0.5 {
        g
    } else {
        1.0 - g
    }
}

pub fn chi2pdf_xq(q: f64, nf: usize) -> f64 {
    chi2_xq(q, nf) / nf as f64
}

pub fn chi2_xq(q: f64, nf: usize) -> f64 {
    let x1 = 0.0;
    let mut x2 = 0.0;
    loop {
        x2 += nf as f64;
        let qt = chi2_df(x2, nf);
        if qt >= q {
            break;
        }
    }
    let f = |x: f64| chi2_df(x, nf);
    fi1(f, q, x1, x2)
}

///C INVERSE OF THE FUNCTION F.
///C RESULT:     FI1=X SUCH THAT Y=F(X).
///C PRECISSION: EPS=1/10**8 can be changed, see below.
///C METHOD:     BISECTING INTERVAL, INITIAL RANGE [X1,X2] (X1<X2).
///C ASSUMPTION: F(X) MONOTON (IN- OR DECREASING) IN THE INTERVAL [X1,X2].
fn fi1<F>(f: F, y: f64, x1: f64, x2: f64) -> f64
where
    F: Fn(f64) -> f64,
{
    let eps = 1.0 / 10f64.powi(8);
    const ITERMAX: usize = 1000;

    let mut x1 = x1;
    let mut x2 = x2;
    let mut y1 = f(x1);
    let mut y2 = f(x2);
    if y1 > y2 {
        (x1, x2, y1, y2) = (x2, x1, y2, y1);
    }

    if y >= y1 && y <= y2 {
        for _ in 0..ITERMAX {
            let x = 0.5 * (x1 + x2);
            if f(x) <= y {
                x1 = x;
            } else {
                x2 = x;
            }
            if (x2 - x1).abs() <= eps {
                return 0.5 * (x2 + x1);
            }
        }
    }

    println!("fi1: no convergence or y out of range !");
    println!("ITERMAX,EPS: {ITERMAX}, {eps}");
    println!("y1,y,y2: {y1}, {y}, {y2}");
    panic!();
}
