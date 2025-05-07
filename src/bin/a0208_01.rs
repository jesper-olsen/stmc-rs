use gnuplot::{AutoOption, AxesCommon, Caption, Color, Figure, PointSymbol};
use stmc_rs::fitl::{LFit, fit_graph, fit_l};

fn main() {
    lfit(true);
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

fn lfit(plot: bool) -> LFit {
    let data: Vec<(f64, f64, f64)> = DATA
        .iter()
        //.map(|(x, y, ey)| subl_1ox(*x, *y, *ey))
        .map(|(x, y, ey)| (*x, *y, *ey))
        .collect();

    let r = fit_l(&data);
    println!("{r}");

    if plot {
        let ([ev1, ev2], [v1, v2, v3, v4]) = fit_graph(&data, &r, -1.0);

        let mut fg = Figure::new();
        fg.set_title("Confidence Ellipse");
        fg.axes2d()
            .set_x_range(AutoOption::Fix(0.65), AutoOption::Fix(1.60))
            .set_y_range(AutoOption::Fix(0.70), AutoOption::Fix(1.35))
            .set_x_label("a_1", &[])
            .set_y_label("a_2", &[])
            .lines(
                ev1.iter().map(|(x, _)| *x).collect::<Vec<f64>>(),
                ev1.iter().map(|(_, y)| *y).collect::<Vec<f64>>(),
                &[Caption("Ellipse"), Color(gnuplot::RGBString("red"))],
            )
            .lines(
                ev2.iter().map(|(x, _)| *x).collect::<Vec<f64>>(),
                ev2.iter().map(|(_, y)| *y).collect::<Vec<f64>>(),
                &[],
            );
        fg.show().unwrap();

        let mut fg = Figure::new();
        fg.set_title("Linear Fit");
        fg.axes2d()
            //.set_x_range(AutoOption::Fix(0.65), AutoOption::Fix(1.60))
            //.set_y_range(AutoOption::Fix(0.70), AutoOption::Fix(1.35))
            .set_x_label("x", &[])
            .set_y_label("y", &[])
            .lines(&v1, &v2, &[Caption(r"fit+e")])
            .lines(&v1, &v3, &[Caption(r"fit")])
            .lines(&v1, &v4, &[Caption(r"fit-e")])
            .y_error_bars(
                DATA.iter().map(|(x, _, _)| *x).collect::<Vec<f64>>(),
                DATA.iter().map(|(_, y, _)| *y).collect::<Vec<f64>>(),
                DATA.iter().map(|(_, _, ey)| *ey).collect::<Vec<f64>>(),
                &[
                    Caption(r"y\_error\_bars"),
                    PointSymbol('T'),
                    Color(gnuplot::RGBString("blue")),
                ],
            );
        fg.show().unwrap();
    }

    r
}

#[cfg(test)]
mod tests {
    use crate::{LFit, lfit};
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

        let r = lfit(false);
        assert_eq!(r, expected);
    }
}
