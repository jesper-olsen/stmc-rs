use crate::gamma::gamma_ln;

pub fn beta(a: f64, b: f64) -> f64 {
    (gamma_ln(a) + gamma_ln(b) - gamma_ln(a + b)).exp()
}

///Incomplete beta function
pub fn beta_i(x: f64, a: f64, b: f64) -> f64 {
    const EPS: f64 = 1.0e-10;
    const ITER_MAX: usize = 200;

    assert!(
        (0.0..=1.0).contains(&x),
        "beta_i: Value is not within the range [0.0, 1.0]"
    );
    let bt = if x == 0.0 || x == 1.0 {
        0.0
    } else {
        (gamma_ln(a + b) - gamma_ln(a) - gamma_ln(b) + a * x.ln() + b * (1.0 - x).ln()).exp()
    };

    let (xx, aa, bb) = if x < (a + 1.0) / (a + b + 2.0) {
        (x, a, b)
    } else {
        (1.0 - x, b, a)
    };
    let apb = aa + bb;
    let ap1 = aa + 1.0;
    let am1 = aa - 1.0;
    let mut bcfm = 1.0;
    let mut bm = 1.0;
    let mut bcf = 1.0;
    let mut bz = 1.0 - apb * xx / ap1;

    for i in 1..=ITER_MAX {
        let xi = i as f64;
        let two_i = xi + xi;
        let c1 = xi * (bb - xi) * xx / ((am1 + two_i) * (aa + two_i));
        let c2 = -(aa + xi) * (apb + xi) * xx / ((aa + two_i) * (ap1 + two_i));
        let bcfp = bcf + c1 * bcfm;
        let bp = bz + c1 * bm;
        let bpp = bp + c2 * bz;
        let bcfold = bcf;
        bcfm = bcfp / bpp;
        bm = bp / bpp;
        bcf = (bcfp + c2 * bcf) / bpp;
        bz = 1.0;

        if (bcf - bcfold).abs() < EPS * bcf.abs() {
            return if x < (a + 1.0) / (a + b + 2.0) {
                bt * bcf / a
            } else {
                1.0 - bt * bcf / b
            };
        }
    }
    panic!("STOP 'BETA_I: A or B too big, or ITER_MAX too small");
}
