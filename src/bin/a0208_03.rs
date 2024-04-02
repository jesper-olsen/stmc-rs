use marsaglia_rs::fitl::{fit_l, subl_power, LFit};

fn main() {
    lfit();
}

// Gaussian data (Bhanot): x, y, error bars
static DATA: [(f64, f64, f64); 5] = [
    (4.0, 0.87739E-01, 0.5E-5),
    (5.0, 0.60978E-01, 0.5E-5),
    (6.0, 0.45411E-01, 0.5E-5),
    (8.0, 0.28596E-01, 0.5E-5),
    (10.0, 0.19996E-01, 0.5E-5),
];

//LINEAR FIT: FIT OF INPUT DATA BY STRAIGHT LINE Y=A*X+B.
fn lfit() -> LFit {
    let data: Vec<(f64, f64, f64)> = DATA
        .iter()
        .map(|(x, y, ey)| subl_power(*x, *y, *ey))
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
            a: (-0.19048445997198254, -1.61852878311417),
            sga: (0.0002805730881271084, 0.00017754271380420746),
            chi2: 1408.4793498125837,
            q: 0.0,
            cov: [
                [0.00000007872125778118214, -0.000000049272544741736555],
                [-0.000000049272544741736555, 0.00000003152141522496272],
            ],
        };

        let r = lfit();
        assert_eq!(r, expected);
    }
}
