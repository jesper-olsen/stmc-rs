use plotters::prelude::*;
use core::ops::Range;

pub fn plot(
    fname: &str,
    caption: &str,
    x_desc: &str,
    y_desc: &str,
    graphs: Vec<(String, Vec<(f64, f64)>)>,
    xrange: Range<f64>,
    yrange: Range<f64>,
) {
    println!("Saving {fname}");
    let root_area = BitMapBackend::new(fname, (600, 400)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(caption, ("sans-serif", 40))
        .build_cartesian_2d(xrange, yrange)
        .unwrap();

    ctx.configure_mesh()
        .x_desc(x_desc)
        .y_desc(y_desc)
        .draw()
        .unwrap();

    graphs.into_iter().enumerate().for_each(|(i, (label, h))| {
        let colour = match i {
            0 => RED,
            1 => BLUE,
            2 => GREEN,
            _ => YELLOW,
        };

        //TODO: generalise LineSeries, AreaSeries
        //.draw_series(LineSeries::new(
        //    x.iter().cloned().zip(uni_cdf.iter().cloned()),
        //    GREEN,
        //))
        ctx.draw_series(
            AreaSeries::new(
                h,
                0.0,             // Baseline
                colour.mix(0.2), // Make the series opac
            )
            .border_style(colour), // Make a brighter border
        )
        .unwrap()
        .label(label)
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], colour));
    });

    ctx.configure_series_labels()
        .border_style(BLACK)
        .background_style(WHITE.mix(0.8))
        .position(SeriesLabelPosition::UpperRight)
        .margin(20)
        .draw()
        .unwrap();
}

pub fn plot2(
    fname: &str,
    caption: &str,
    x_desc: &str,
    y_desc: &str,
    graphs: Vec<(String, Vec<(f64, f64)>)>,
    xrange: Range<f64>,
    yrange: Range<f64>,
) {
    println!("Saving {fname}");
    let root_area = BitMapBackend::new(fname, (900, 600)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(caption, ("Arial", 40))
        .build_cartesian_2d(xrange, yrange)
        .unwrap();

    ctx.configure_mesh()
        .x_desc(x_desc)
        .y_desc(y_desc)
        .draw()
        .unwrap();

    graphs.into_iter().for_each(|(_, h)| {
        ctx.draw_series(LineSeries::new(h, BLACK)).unwrap();
    });
}
