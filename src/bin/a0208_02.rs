use marsaglia_rs::fitl::{fit_l, subl_1ox, LFit};

fn main() {
    lfit();
}

// Gaussian data: x, y, error bars
static DATA: [(f64, f64, f64); 6] = [
    (16.0, 0.1086, 0.0007),
    (24.0, 0.1058, 0.0008),
    (34.0, 0.1039, 0.0012),
    (50.0, 0.1016, 0.0008),
    (70.0, 0.0995, 0.0012),
    (100.0, 0.0990, 0.0012),
];

//LINEAR FIT: FIT OF INPUT DATA BY STRAIGHT LINE Y=A*X+B.
fn lfit() -> LFit {
    let data: Vec<(f64, f64, f64)> = DATA
        .iter()
        .map(|(x, y, ey)| subl_1ox(*x, *y, *ey))
        .collect();

    let r = fit_l(&data);
    println!("{r}");

    r
}

#[cfg(test)]
mod tests {
    use crate::{lfit, LFit};
    #[test]
    fn lfit_test() {
        let expected = LFit {
            a: (0.17829769783513708, 0.09782177231480599),
            sga: (0.0192082135777534, 0.0007854559197949675),
            chi2: 2.07554626333672,
            q: 0.7218661665309511,
            cov: [
                [0.00036895546884859, -0.000013293884969307767],
                [-0.000013293884969307767, 0.0000006169410019409585],
            ],
        };

        let r = lfit();
        assert_eq!(r, expected);
    }
}
