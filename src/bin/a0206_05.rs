use marsaglia_rs::kolm::{kolm1, kolm2_as};
use marsaglia_rs::marsaglia::Marsaglia;

/// Comparision of the asymptotic 2-sided Kolmogorov test with
/// results from the exact one sided test.

fn kol2() -> Vec<(f64, f64, f64)> {
    let mut l = Vec::new();
    const N: usize = 5;
    const NITER: usize = 4;
    let mut data = [0.0f64; N];
    let mut fxct = [0.0f64; N];

    let mut rng = Marsaglia::new(12, 34, 56, 78);
    println!(" KOLM2_AS with N={N}");
    for iter in 0..NITER {
        for e in data.iter_mut() {
            *e = rng.uni();
        }
        data.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (e1, e2) in fxct.iter_mut().zip(data.iter()) {
            *e1 = *e2;
        }
        let (_del1, _del2, q1, q2) = kolm1(&fxct);
        let q3 = q1.min(q2);
        let (del, q) = kolm2_as(&fxct);
        println!(
            "ITER,del,q3,2q3,q_as: {iter:4} {del:10.5} {q3:10.5} {twoq3:10.5} {q:10.5}",
            twoq3 = 2.0 * q3
        );
        l.push((del, q3, q));
    }
    l
}

fn main() {
    kol2();
}

#[cfg(test)]
mod tests {
    use crate::kol2;

    #[test]
    fn kol1_test() {
        let r: Vec<(f64, f64, f64)> = vec![
            (0.28297039270401003, 0.3824337983202357, 0.7432117198891041),
            (0.3501267910003662, 0.2364588367221349, 0.4773271312715121),
            (0.32357357740402226, 0.29054150885185837, 0.5798517059525516),
            (
                0.5319563031196595,
                0.037752917264718334,
                0.07568942457671027,
            ),
        ];
        let r2 = kol2();
        assert_eq!(r, r2);
    }
}
