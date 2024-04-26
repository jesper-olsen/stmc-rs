use stmc_rs::cau::cau_df;
use stmc_rs::kolm::kolm1;
use stmc_rs::marsaglia::Marsaglia;

// Examples for the 1-sided Kolmogorov test.

fn kol1() -> Vec<(usize, f64, f64)> {
    let mut l = Vec::new();
    for lcau in [true, false] {
        let mut data = [0.0; 100];
        let mut fxct = [0.0; 100];
        println!(
            "{} random numbers:",
            if lcau { "Cauchy" } else { "Uniform" }
        );
        for n in [8, 100] {
            let mut rng = Marsaglia::new(12, 34, 56, 78);
            for e in data.iter_mut().take(n) {
                *e = if lcau { rng.cauchy() } else { rng.uni() };
            }
            data[..n].sort_by(|a, b| a.partial_cmp(b).unwrap());
            for i in 0..n {
                fxct[i] = if lcau { cau_df(data[i]) } else { data[i] };
            }
            let (_del1, _del2, q1, q2) = kolm1(&fxct[..n]);
            println!("N: {n:6} Q1: {q1:10.2} Q2: {q2:10.2}");
            l.push((n, q1, q2));
        }
    }
    l
}

fn main() {
    kol1();
}

#[cfg(test)]
mod tests {
    use crate::kol1;

    #[test]
    fn kol1_test() {
        let r: Vec<(usize, f64, f64)> = vec![
            (8, 0.2736615841466952, 0.6609602759550518),
            (100, 0.2207449986218413, 0.23087660796738066),
            (8, 0.2938562438611388, 0.9340610497326689),
            (100, 0.13920192741340637, 0.7980354361862717),
        ];

        let r2 = kol1();
        assert_eq!(r, r2);
    }
}
