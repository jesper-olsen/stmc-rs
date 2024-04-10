use crate::gamma::gamma_p;
use std::f64::consts::PI;
use std::fmt;

// LINEAR FIT: FIT OF INPUT DATA BY STRAIGHT LINE Y=A*X+B.
//
// INPUT:   GAUSSIAN DATA          Y (I), I=1,...,N
//          AND THEIR ERROR BARS   YE(I), I=1,...,N.
//
// OUTPUT:  CONSTANTS  A(1), A(2) WITH ERRORS BARS SGA(1), SGA(2).
//          CORRESPONDING TO THE FIT Y=A(1)+A(2)*X. IF N.GE.3:
//          ALSO  CHI2  AND  THE GOODNESS OF FIT Q ARE RETURNED.

#[derive(Debug, PartialEq)]
pub struct LFit {
    pub a: (f64, f64),
    pub sga: (f64, f64),
    pub chi2: f64,
    pub q: f64,
    pub cov: [[f64; 2]; 2],
}

impl fmt::Display for LFit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "LFit:")?;
        writeln!(f, "A(1) = {} +/- {}", self.a.0, self.sga.0)?;
        writeln!(f, "A(2) = {} +/- {}", self.a.1, self.sga.1)?;
        writeln!(f, "\nCOVARIANCE MATRIX:")?;
        writeln!(f, "{} {}", self.cov[0][0], self.cov[0][1])?;
        writeln!(f, "{} {}", self.cov[1][0], self.cov[1][1])?;
        writeln!(f, "\nSTATISTICAL ANALYSIS:")?;
        writeln!(f, "CHI2 = {}", self.chi2)?;
        writeln!(f, "   Q = {}", self.q)
    }
}

pub fn fit_l(data: &[(f64, f64, f64)]) -> LFit {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut s = 0.0;
    for (x, y, sgy) in data {
        let wt = 1.0 / sgy.powi(2);
        s += wt;
        sx += x * wt;
        sy += y * wt;
    }
    let sxos = sx / s;

    let mut stt = 0.0;
    let mut a = (0.0, 0.0);
    for (x, y, sgy) in data {
        let t = (x - sxos) / sgy;
        stt += t * t;
        a.1 += t * y / sgy;
    }
    a.1 /= stt;
    a.0 = (sy - sx * a.1) / s;

    let chi2 = data
        .iter()
        .map(|(x, y, sgy)| ((y - a.0 - a.1 * x) / sgy).powi(2))
        .sum();
    let q = if data.len() > 2 {
        1.0 - gamma_p(0.5 * (data.len() - 2) as f64, 0.5 * chi2)
    } else {
        0.0
    };
    let sga = (
        ((1.0 + sx.powi(2) / (s * stt)) / s).sqrt(),
        (1.0 / stt).sqrt(),
    );
    let cov = -sx / (s * stt);

    LFit {
        a,
        sga,
        chi2,
        q,
        cov: [[sga.0.powi(2), cov], [cov, sga.1.powi(2)]],
    }
}

// Fit y=c1*exp(-c2*ln(x)) -> ynew=ln(y)=ln(c1)-c2*xnew
// with xnew=ln(x). Then: a1=ln(c1), a2=-c2.
pub fn subl_power(x: f64, y: f64, ey: f64) -> (f64, f64, f64) {
    let yup = (y + ey).ln();
    let yt = (yup + (y - ey).ln()) / 2.0;
    (x.ln(), yt, yup - yt)
}

// Fit y=c1+c2/x -> ynew=c1*x+c2=a1+a2*x
pub fn subl_1ox(x: f64, y: f64, ey: f64) -> (f64, f64, f64) {
    (x, y * x, ey * x)
}

const M2X2_1: [[f64; 2]; 2] = [[1.0, 0.0], [0.0, 1.0]];
const EPS: f64 = 1e-10;

pub fn eigen_2x2(amat: &[[f64; 2]; 2]) -> ([f64; 2], [[f64; 2]; 2]) {
    let phalf = -0.5 * (amat[0][0] + amat[1][1]);
    let q = amat[0][0] * amat[1][1] - amat[0][1] * amat[1][0];
    let arg = phalf.powi(2) - q;
    if arg < 0.0 {
        panic!("eigen_2x2 cannot be used; Eigenvalues are complex!");
    }
    let eval = [-phalf + arg.sqrt(), -phalf - arg.sqrt()];
    let mut evct = [[0.0f64; 2]; 2];
    for i in 0..=1 {
        let deno0 = eval[i] - amat[0][0];
        let deno1 = eval[i] - amat[1][1];
        if deno0.abs() > deno1.abs() {
            if deno0.abs() == EPS {
                return (eval, M2X2_1);
            } else {
                let factor = amat[0][1] / deno0;
                evct[1][i] = 1.0 / (1.0 + factor.powi(2)).sqrt();
                evct[0][i] = factor * evct[1][i];
            }
        } else if deno1.abs() == EPS {
            return (eval, M2X2_1);
        } else {
            let factor = amat[1][0] / deno1;
            evct[0][i] = 1.0 / (1.0 + factor.powi(2)).sqrt();
            evct[1][i] = factor * evct[0][i];
        }
    }
    (eval, evct)
}

// draw confidence ellipse
fn ellipse(
    xm: (f64, f64),
    sgxm: (f64, f64),
    sgym: (f64, f64),
    eigvs: &[[f64; 2]; 2],
    prob: f64,
) -> [Vec<(f64, f64)>; 2] {
    const NPOINTS: usize = 100;
    // 1. ELLIPSE:
    // C PROB<0: STANDARD COVARIANCE ELLIPSE WITH 39% CONFIDENCE.
    let r = if prob < 0.0 {
        1.0
    } else {
        (-2.0 * (1.0 - prob).ln()).sqrt()
    };
    let delphy = 2.0 * PI / NPOINTS as f64;
    let mut v1 = Vec::new();
    for i in 0..=NPOINTS {
        let phy = i as f64 * delphy;
        let y1 = r * sgym.0 * phy.sin();
        let y2 = r * sgym.1 * phy.cos();
        let x1 = xm.0 + eigvs[0][0] * y1 + eigvs[0][1] * y2;
        let x2 = xm.1 + eigvs[1][0] * y1 + eigvs[1][1] * y2;
        v1.push((x1, x2));
    }

    //C 2. STANDARD DEVIATIONS IN ORIGNAL DIRECTIONS:
    let v2 = vec![
        (xm.0 + sgxm.0, xm.1),
        (xm.0 - sgxm.0, xm.1),
        (xm.0, xm.1),
        (xm.0, xm.1 + sgxm.1),
        (xm.0, xm.1 - sgxm.1),
        (xm.0, xm.1),
        //C 3. STANDARD DEVIATIONS IN THE INDEPENDENT DIRECTIONS:
        //println!("SET TEXTURE DASHES");
        (xm.0 + sgym.0 * eigvs[0][0], xm.1 + sgym.0 * eigvs[1][0]),
        (xm.0 - sgym.0 * eigvs[0][0], xm.1 - sgym.0 * eigvs[1][0]),
        (xm.0, xm.1),
        (xm.0 + sgym.1 * eigvs[0][1], xm.1 + sgym.1 * eigvs[1][1]),
        (xm.0 - sgym.1 * eigvs[0][1], xm.1 - sgym.1 * eigvs[1][1]),
        (xm.0, xm.1),
    ];

    [v1, v2]
}

pub fn fit_graph(
    data: &[(f64, f64, f64)],
    lfit: &LFit,
    prob: f64,
) -> ([Vec<(f64, f64)>; 2], [Vec<f64>; 4]) {
    const NFIG: usize = 200;

    let (covdd, eigvs) = eigen_2x2(&lfit.cov);
    let sgb = (covdd[0].sqrt(), covdd[1].sqrt());
    let egraph = ellipse(lfit.a, lfit.sga, sgb, &eigvs, prob);

    // 3. INCLUDE CONFIDENCE LIMITS FOR REGRESSION LINE:
    let (xn, _, _) = data.last().unwrap();
    let (x0, _, _) = data.first().unwrap();
    let delx = (xn - x0) / NFIG as f64;

    let mut v1 = Vec::new();
    let mut v2 = Vec::new();
    let mut v3 = Vec::new();
    let mut v4 = Vec::new();
    for i in 0..=NFIG {
        let xx = x0 + i as f64 * delx;
        let varb1 = (eigvs[0][0] + xx * eigvs[1][0]).abs().powi(2) * covdd[0];
        let varb2 = (eigvs[0][1] + xx * eigvs[1][1]).abs().powi(2) * covdd[1];
        let ym = lfit.a.0 + lfit.a.1 * xx - (varb1 + varb2).sqrt();
        let yy = lfit.a.0 + lfit.a.1 * xx;
        let yp = lfit.a.0 + lfit.a.1 * xx + (varb1 + varb2).sqrt();
        v1.push(xx);
        v2.push(ym);
        v3.push(yy);
        v4.push(yp);
    }

    (egraph, [v1, v2, v3, v4])
}
