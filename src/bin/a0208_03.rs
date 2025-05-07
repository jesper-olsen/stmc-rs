use gnuplot::{AxesCommon, Caption, Color, Figure, PointSymbol};
use stmc_rs::fitl::{LFit, fit_graph, fit_l, subl_power, subplot_power};

fn main() {
    lfit(true);
}

// bhanot.D - x, y, error bars
static DATA: [(f64, f64, f64); 5] = [
    (4.0, 0.87739E-01, 0.5E-05),
    (5.0, 0.60978E-01, 0.5E-05),
    (6.0, 0.45411E-01, 0.5E-05),
    (8.0, 0.28596E-01, 0.5E-05),
    (10.0, 0.19996E-01, 0.5E-05),
];

//LINEAR FIT: FIT OF INPUT DATA BY STRAIGHT LINE Y=A*X+B.
fn lfit(plot: bool) -> LFit {
    let data: Vec<(f64, f64, f64)> = DATA
        .iter()
        .map(|(x, y, ey)| subl_power(*x, *y, *ey))
        .collect();

    let r = fit_l(&data);
    println!("{r}");

    if plot {
        let (vx, vy) = subplot_power(&r, &data, 200);
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
                &[
                    Caption(r"y\_error\_bars"),
                    PointSymbol('T'),
                    Color(gnuplot::RGBString("blue")),
                ],
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
            a: (-0.19048445997198254, -1.61852878311417),
            sga: (0.0002805730881271084, 0.00017754271380420746),
            chi2: 1408.4793498125837,
            q: 0.0,
            cov: [
                [0.00000007872125778118214, -0.000000049272544741736555],
                [-0.000000049272544741736555, 0.00000003152141522496272],
            ],
        };

        let r = lfit(false);
        assert_eq!(r, expected);
    }
}
