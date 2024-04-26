use gnuplot::{AxesCommon, Caption, Color, Figure, PointSymbol};
use stmc_rs::fitl::{fit_graph, fit_l, LFit, subl_1ox, subplot_1ox};

fn main() {
    lfit(true);
}

// BERG.D - x, y, error bars
static DATA: [(f64, f64, f64); 6] = [
    (16.0, 0.1086, 0.0007),
    (24.0, 0.1058, 0.0008),
    (34.0, 0.1039, 0.0012),
    (50.0, 0.1016, 0.0008),
    (70.0, 0.0995, 0.0012),
    (100.0, 0.0990, 0.0012),
];

fn lfit(plot: bool) -> LFit {
    let data: Vec<(f64, f64, f64)> = DATA
        .iter()
        .map(|(x, y, ey)| subl_1ox(*x, *y, *ey))
        //.map(|(x, y, ey)| (x, y, ey))
        .collect();

    let r = fit_l(&data);
    println!("{r}");

    if plot {
        let (vx, vy) = subplot_1ox(&r, &data, 200);
        let mut fg = Figure::new();
        fg.set_title("Linear Fit");
        fg.axes2d()
            .set_x_label("x", &[])
            .set_y_label("y", &[])
            .lines(&vx, &vy, &[Caption(r"fit")])
            .y_error_bars(
                DATA.iter().map(|(x, _, _)| *x).collect::<Vec<f64>>(),
                DATA.iter().map(|(_, y, _)| *y).collect::<Vec<f64>>(),
                DATA.iter().map(|(_, _, ey)| *ey).collect::<Vec<f64>>(),
                &[Caption(r"y\_error\_bars"), PointSymbol('T'), Color("blue")],
            );
        fg.show().unwrap();

        let ([ev1, ev2], [v1, v2, v3, v4]) = fit_graph(&data, &r, -1.0);

        let mut fg = Figure::new();
        fg.set_title("Confidence Ellipse");
        fg.axes2d()
            //.set_x_range(AutoOption::Fix(0.65), AutoOption::Fix(1.60))
            //.set_y_range(AutoOption::Fix(0.70), AutoOption::Fix(1.35))
            .set_x_label("a_1", &[])
            .set_y_label("a_2", &[])
            .lines(
                ev1.iter().map(|(x, _)| *x).collect::<Vec<f64>>(),
                ev1.iter().map(|(_, y)| *y).collect::<Vec<f64>>(),
                &[Caption("Ellipse"), Color("red")],
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
            .set_x_label("x", &[])
            .set_y_label("y", &[])
            .lines(&v1, &v2, &[Caption(r"fit+e")])
            .lines(&v1, &v3, &[Caption(r"fit")])
            .lines(&v1, &v4, &[Caption(r"fit-e")])
            .y_error_bars(
                DATA.iter().map(|(x, _, _)| *x).collect::<Vec<f64>>(),
                DATA.iter().map(|(_, y, _)| *y).collect::<Vec<f64>>(),
                DATA.iter().map(|(_, _, ey)| *ey).collect::<Vec<f64>>(),
                &[Caption(r"y\_error\_bars"), PointSymbol('T'), Color("blue")],
            );
        fg.show().unwrap();
    }

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

        let r = lfit(false);
        assert_eq!(r, expected);
    }
}
