use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::plot::plot2;
use marsaglia_rs::steb::{bias, datjack, steb0, stebj0};

// Copyright, Bernd Berg, Oct 31, 2000.
// VARIANCE EXAMPLE FOR DOUBLE JACKKNIFE BIAS CORRECTED
// ESTIMATORS (Berg, Comp. Phys. 69 (92) 7-14).

fn main() {
    bias_var();
}

const N: usize = 16;
const N1: usize = N - 1;

fn bias_var() -> (f64,f64) {
    let mut data = [0.0f64; N];
    let mut dat2 = [0.0f64; N];
    let mut datj = [0.0f64; N];
    let mut datj2 = [0.0f64; N];
    let mut datjj = [[0.0f64; N1]; N];
    let mut datjj2 = [[0.0f64; N1]; N];
    let mut datj2mm = [0.0f64; N];

    let mut rng = Marsaglia::new(12, 34, 56, 78);
    for (e1, e2) in data.iter_mut().zip(dat2.iter_mut()) {
        *e1 = rng.uni();
        *e2 = (*e1).powi(2);
    }
    println!("NUMERICAL ESTIMATES FOR THE MEAN SQUARED:\n");
    let datm = data.iter().sum::<f64>() / data.len() as f64;
    let datm2 = datm.powi(2);
    println!(" BIASED ESTIMATE FROM ALL DATA: {datm2}");

    let (badm2, badm2e) = steb0(&dat2);
    println!(" BAD, BIASED ESTIMATOR (FALSE): {badm2} +/- {badm2e}");

    let da2m = dat2.iter().sum::<f64>() / dat2.len() as f64;
    let datmm2 = datm2 - (da2m - datm2) / (dat2.len() - 1) as f64;
    println!("UNBIASED ESTIMATOR:             {datmm2}");

    datjack(&data, &mut datj);
    for (e1, e2) in datj2.iter_mut().zip(datj.iter()) {
        *e1 = (*e2).powi(2);
    }
    let (datj2m, _, error0) = stebj0(&datj2);
    println!(" STANDARD, BIASED JACKNIFE ESTIMATOR: {datj2m} +/- {error0}");

    datjack2(&data, &mut datjj);
    for i1 in 0..N {
        for i2 in 0..N1 {
            datjj2[i1][i2] = datjj[i1][i2].powi(2);
        }
    }
    let (datmm2, _, error) = stebjj1(&datjj2, &datj2, &mut datj2mm);

    println!(" SECOND LEVEL");
    println!(" BIAS-CORRECTED JACKKNIFE ESTIMATOR: {datmm2} +/- {error}");
    (datmm2,error)
}

// CALCULATION OF  SECOND LEVEL JACKKNIFE BINS  XJJ(N-1,N)
// FROM  N  DATA  X(N).
fn datjack2(x: &[f64; N], xjj: &mut [[f64; N1]; N]) {
    let xsum = x.iter().sum::<f64>();

    let nm1 = N - 1;
    let factor = 1.0 / (N - 2) as f64;

    for i in 0..N {
        for j in 0..nm1 {
            let mut jj = j;
            if j >= i {
                jj += 1;
            }
            xjj[i][j] = factor * (xsum - x[i] - x[jj]);
        }
    }
}

// CALCULATION OF JACKKNIFE ESTIMATORS FROM INPUT OF SECOND
// FJJ(N-1,N)  AND FIRST  FJ(N)  LEVEL JACKKNIFE BINS:
//
// - THE JACKKNIFE SAMPLE OF BIAS CORRECTED MEAN VALUES FJMM(N),
// - THE BIAS CORRECTED MEAN VALUE FMM,
// - THE VARIANCE  FV  AND THE ERROR BAR  FE  FOR  THE
//   BIAS CORRECTED MEAN VALUE.
pub fn stebjj1(fjj: &[[f64; N1]; N], fj: &[f64; N], fjmm: &mut [f64; N]) -> (f64, f64, f64) {
    assert!(fj.len() > 2, "stebjj1: too small");
    assert!(fjj.len() == fj.len(), "stebjj1: mismatched dimensions");
    for i in 0..fj.len() {
        (fjmm[i], _) = bias(&fjj[i], fj[i]);
    }
    stebj0(fjmm)
}

#[cfg(test)]
mod tests {
    use crate::bias_var;
    #[test]
    fn bias_var() {
        let tup = bias_var();
        assert_eq!(bias, (0.23837126846653026 ,0.07978421694955856));
    }
}
