use marsaglia_rs::fitl::{fit_l, subl_1ox, LFit, eigen_2x2};
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

    fit_lgnu("fig_a0208_01", &data, &r, -1.0);
    r
}

// draw confidence ellipse
fn ellipse(fname: &str, xm: (f64, f64), sgxm: (f64, f64), sgym: (f64, f64), eigvs: &[[f64; 2]; 2], prob: f64) {
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
    let s=format!("{fname}a.png");       
    plot2(
        s.as_str(),
        "Confidence Ellipse", // title
        "a_1",
        "a_2",
        graphs,
        0.65..1.6,
        0.7..1.35,
    );
}

fn fit_lgnu(fname: &str, data: &[(f64, f64, f64)], lfit: &LFit, prob: f64) {
    const NFIG: usize = 200;

    // lfit.d1
    //for (x,y,ey) in data {
    //    println!("{x} {y} {ey}");
    //}

    let (covdd, eigvs) = eigen_2x2(&lfit.cov);
    let sgb = (covdd[0].sqrt(), covdd[1].sqrt());
    ellipse(fname, lfit.a, lfit.sga, sgb, &eigvs, prob);

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
    let s=format!("{fname}b.png");       
    plot2(
        s.as_str(),
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
