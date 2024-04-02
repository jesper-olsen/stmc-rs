use marsaglia_rs::fitl::{fit_graph, fit_l, LFit};
use marsaglia_rs::plot::plot2;

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

    let [egraph, lgraph] = fit_graph(&data, &r, -1.0);

    plot2(
        "fig_a0208_01a.png",
        "Confidence Ellipse", // title
        "a_1",
        "a_2",
        egraph,
        0.65..1.6,
        0.7..1.35,
    );
    plot2(
        "fig_a0208_01b.png",
        "lfit", // title
        "",
        "",
        lgraph,
        0.0..3.0,
        0.0..5.0,
    );

    r
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
