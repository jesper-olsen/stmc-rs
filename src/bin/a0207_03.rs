use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::steb::{bias, datjack, steb0, stebj0};

fn main() {
    bias_big();
}

const N: usize = 32;
const N1: usize = N - 1;

// ILLUSTRATION OF THE DOUBLE JACKKNIFE BIAS CORRECTED ESTIMATOR.
// Here the bias is big through the ABS in the EXP function.
// (Related paper: Berg, Comp. Phys. Commun. 69 (92) 7-14.)
fn bias_big() -> (f64, f64) {
    //const NTRY: usize = 2;
    const NTRY: usize = 2000;
    const NPR: usize = 4;
    const A0: f64 = 2.0;
    const NDAT: usize = 320;
    let mut stup2m = [0.0f64; NTRY];
    let mut datm2 = [0.0f64; NTRY];
    let mut datj2m = [0.0f64; NTRY];
    let mut datmm2 = [0.0f64; NTRY];
    let mut dat0 = [0.0f64; NDAT];
    let mut data = [0.0f64; N];
    let mut stup2 = [0.0f64; N];
    let mut dat2 = [0.0f64; N];
    let mut datj = [0.0f64; N];
    let mut datj2 = [0.0f64; N];
    let mut datjj = [[0.0f64; N1]; N];
    let mut datjj2 = [[0.0f64; N1]; N];
    let mut datj2mm = [0.0f64; N];

    println!("EXACT: {}", 1.0 / 3.0);
    let mut rng = Marsaglia::new(12, 34, 56, 78);

    for itry in 0..NTRY {
        if itry < NPR {
            println!("\nITRY = {}", itry + 1);
        }
        for e in dat0.iter_mut() {
            *e = rng.gauss();
            _ = rng.gauss();
        }
        bining(&dat0, &mut data);
        for e in dat0.iter_mut() {
            *e = (*e).powi(4);
        }
        bining(&dat0, &mut dat2);
        let datm = data.iter().sum::<f64>() / data.len() as f64;
        let da2m = dat2.iter().sum::<f64>() / dat2.len() as f64;
        datm2[itry] = (A0 * datm.abs()).exp() / da2m;
        if itry < NPR {
            println!(" BAD, BIASED ESTIMATOR FROM ALL DATA: {}", datm2[itry]);
        }

        datjack(&data, &mut datj);
        datjack(&dat2, &mut datj2);
        for i in 0..N {
            stup2[i] = (A0 * data[i].abs()).exp() / dat2[i];
            datj2[i] = (A0 * datj[i].abs()).exp() / datj2[i];
        }

        let error00;
        (stup2m[itry], error00) = steb0(&stup2);
        let error01;
        (datj2m[itry], _, error01) = stebj0(&datj2);
        if itry <= NPR {
            println!(" BAD, BIASED ESTIMATOR: {} +/- {error00}", stup2m[itry]);
            println!(
                " STANDARD, BIASED JACKKNIFE ESTIMATOR: {} +/- {error01}",
                datj2m[itry]
            );
        }

        datjack2(&data, &mut datjj);
        datjack2(&dat2, &mut datjj2);
        for i1 in 0..N {
            for i2 in 0..N1 {
                datjj2[i1][i2] = (A0 * datjj[i1][i2].abs()).exp() / datjj2[i1][i2];
            }
        }
        for i in 0..datj2.len() {
            (datj2mm[i], _) = bias(&datjj2[i], datj2[i]);
        }
        let error;
        (datmm2[itry], _, error) = stebj0(&datj2mm);
        if itry < NPR {
            println!(
                "SECOND LEVEL BIAS-CORRECTED JACKKNIFE ESTIMATOR: {} +/- {error}",
                datmm2[itry]
            );
        }
    }

    println!("Average over NTRY = {NTRY} calculations.");
    println!("Comparison of results:");

    let (xm2, em2) = steb0(&datm2);
    println!("BIASED ESTIMATE FROM ALL DATA: {xm2} +/- {em2}:");
    let (xs2m, est2m) = steb0(&stup2m);
    println!("BAD BIASED ESTIMATOR: {xs2m} +/- {est2m}:");

    let (xj2m, ej2m) = steb0(&datj2m);
    println!("STANDARD, BIASED JACKNIFE ESTIMATOR: {xj2m} +/- {ej2m}:");

    let (xmm2, emm2) = steb0(&datmm2);
    println!("SECOND LEVEL BIAS_CORRECTED JACKNIFE ESTIMATOR: {xmm2} +/- {emm2}:");
    (xmm2, emm2)
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

// CALCULATION OF BINNED DATA (datb) FROM  ORIGINAL DATA.
fn bining(data: &[f64], bdata: &mut [f64]) {
    // split into segments (NBINS) and average each segment
    let bin_len = data.len() / bdata.len();
    assert_eq!(bin_len * bdata.len(), data.len());
    for i in 0..bdata.len() {
        bdata[i] = 0.0;
        for j in 0..bin_len {
            bdata[i] += data[i * bin_len + j];
        }
        bdata[i] /= bin_len as f64;
    }
}

#[cfg(test)]
mod tests {
    use crate::bias_big;
    #[test]
    fn test_bias_big() {
        assert_eq!(bias_big(), (0.34666991387440965, 0.0017375655697640011));
    }
}
