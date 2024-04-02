use marsaglia_rs::fitl::{fit_l, subl_1ox, LFit};
use marsaglia_rs::plot::{plot2};
use std::f64::consts::PI;

fn main() {
    lfit();
}

// Gaussian data: x, y, error bars
// Linear Regression. Example 9-2 (p.201) from S. Brandt: Statistical and Computational
// Methods in Data Analysis.
static DATA: [(f64, f64, f64); 4] = [
    (0.0, 1.4, 0.5),
    (1.0, 1.5, 0.5),
    (2.0, 3.7, 0.5),
    (3.0, 4.1, 0.5),
];

fn lfit() -> LFit {
    let data: Vec<(f64, f64, f64)> = DATA
        .iter()
        //.map(|(x, y, ey)| subl_1ox(*x, *y, *ey))
        .map(|(x, y, ey)| (*x, *y, *ey))
        .collect();

    let r = fit_l(&data);
    println!("{r}");

    fit_lgnu(&data, &r, -1.0);
    r
}

const M2X2_1: [[f64; 2]; 2] = [[1.0, 0.0], [0.0, 1.0]];
const EPS: f64 = 1e-10;

fn eigen_2x2(amat: &[[f64; 2]; 2]) -> ([f64; 2], [[f64; 2]; 2]) {
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
        } else {
            if deno1.abs() == EPS {
                return (eval, M2X2_1);
            } else {
                let factor = amat[1][0] / deno1;
                evct[0][i] = 1.0 / (1.0 + factor.powi(2)).sqrt();
                evct[1][i] = factor * evct[0][i];
            }
        }
    }
    (eval, evct)
}

// draw confidence ellipse
fn ellipse(xm: (f64, f64), sgxm: (f64, f64), sgym: (f64, f64), eigvs: &[[f64; 2]; 2], prob: f64) {
    const NGNU: usize = 100;
    // 1. ELLIPSE:
    // C PROB<0: STANDARD COVARIANCE ELLIPSE WITH 39% CONFIDENCE.
    // create ellipse.d1
    println!("ellipse.d1");
    let r = if prob < 0.0 {
        1.0
    } else {
        (-2.0 * (1.0 - prob).ln()).sqrt()
    };
    let delphy = 2.0 * PI / NGNU as f64;
    let mut v1 = Vec::new();
    for i in 0..=NGNU {
        let phy = i as f64 * delphy;
        let y1 = r * sgym.0 * phy.sin();
        let y2 = r * sgym.1 * phy.cos();
        let x1 = xm.0 + eigvs[0][0] * y1 + eigvs[0][1] * y2;
        let x2 = xm.1 + eigvs[1][0] * y1 + eigvs[1][1] * y2;
        v1.push((x1, x2));
    }

    //C 2. STANDARD DEVIATIONS IN ORIGNAL DIRECTIONS:
    // create ellipse.d2
    println!("ellipse.d2");
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

    let graphs = vec![(String::new(), v1), (String::new(), v2)];
    plot2(
        "fig_a0208_01a.png",
        "Confidence Ellipse", // title
        "a_1",
        "a_2",
        graphs,
        0.65..1.6,
        0.7..1.35,
    );
}

fn fit_lgnu(data: &[(f64, f64, f64)], lfit: &LFit, prob: f64) {
    const NFIG: usize = 200;

    // lfit.d1
    //for (x,y,ey) in data {
    //    println!("{x} {y} {ey}");
    //}

    let (covdd, eigvs) = eigen_2x2(&lfit.cov);
    let sgb = (covdd[0].sqrt(), covdd[1].sqrt());
    ellipse(lfit.a, lfit.sga, sgb, &eigvs, prob);

    // 3. INCLUDE CONFIDENCE LIMITS FOR REGRESSION LINE:
    println!("lfit.d2");
    let (xn, _, _) = data.last().unwrap();
    let (x0, _, _) = data.first().unwrap();
    let delx = (xn - x0) / NFIG as f64;
    let mut v1=Vec::new();
    let mut v2=Vec::new();
    let mut v3=Vec::new();
    for i in 0..=NFIG {
        let xx = x0 + i as f64 * delx;
        let varb1 = (eigvs[0][0] + xx * eigvs[1][0]).abs().powi(2) * covdd[0];
        let varb2 = (eigvs[0][1] + xx * eigvs[1][1]).abs().powi(2) * covdd[1];
        let ym = lfit.a.0 + lfit.a.1 * xx - (varb1 + varb2).sqrt();
        let yy = lfit.a.0 + lfit.a.1 * xx;
        let yp = lfit.a.0 + lfit.a.1 * xx + (varb1 + varb2).sqrt();
        //println!("{xx} {ym} {yy} {yp}");
        v1.push((xx,ym));
        v2.push((xx,yy));
        v3.push((xx,yp))
    }
    let graphs = vec![(String::new(), v1), (String::new(), v2), (String::new(),v3)];
    plot2(
        "fig_a0208_01b.png",
        "lfit", // title
        "",
        "",
        graphs,
        0.0..3.0,
        0.0..5.0,
    );
}

#[cfg(test)]
mod tests {
    use crate::{lfit, LFit};
    #[test]
    fn lfit_test() {
        let expected = LFit {
            a: (1.13, 1.03),
            sga: (0.4183300132670378, 0.22360679774997896),
            chi2: 3.1320000000000006,
            q: 0.20887902975523354,
            cov: [
                [0.17500000000000002, -0.075],
                [-0.075, 0.049999999999999996],
            ],
        };

        let r = lfit();
        assert_eq!(r, expected);
    }
}
