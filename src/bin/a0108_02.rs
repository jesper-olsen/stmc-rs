use marsaglia_rs::marsaglia::Marsaglia;
use plotters::prelude::*;
use std::f64::consts::PI;

fn plot(fname: &str, graphs: Vec<(String, Vec<(f64, f64)>)>, xmin: f64, xmax: f64, ymax: f64) {
    println!("Saving {fname}");
    let root_area = BitMapBackend::new(fname, (600, 400)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption("Probability Densities", ("sans-serif", 40))
        .build_cartesian_2d(xmin..xmax, 0f64..ymax)
        .unwrap();

    ctx.configure_mesh().x_desc("x").y_desc("f").draw().unwrap();

    graphs.into_iter().enumerate().for_each(|(i, (label, h))| {
        let colour = match i {
            0 => RED,
            _ => BLUE,
        };
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

fn gaussian_pdf(x: f64, mean: f64, std_dev: f64) -> f64 {
    let exponent = -((x - mean) * (x - mean)) / (2.0 * std_dev * std_dev);
    let coefficient = 1.0 / (std_dev * (2.0 * PI).sqrt());
    coefficient * exponent.exp()
}

fn bin_it(data: &[f64], nbins: usize, bdata: &mut [f64]) {
    // split into segments (NBINS) and average each segment
    let bin_len = data.len() / nbins;
    assert_eq!(bin_len * nbins, data.len());
    for i in 0..nbins {
        for j in 0..bin_len {
            bdata[i] += data[i * bin_len + j];
        }
        bdata[i] /= bin_len as f64;
    }
}

fn main() {
    const NDAT: usize = 10000;
    const NBINS: usize = 500;
    let bin_len = NDAT / NBINS;

    let sdv = 1.0 / (12.0 * bin_len as f64).sqrt();
    let xmin = 0.5 - 4.0 * sdv;
    let xmax = 0.5 + 4.0 * sdv;
    println!("xmin: {xmin} xmax: {xmax}");

    let mut rng = Marsaglia::new(12, 34, 56, 78);
    let mut data = [0.0f64; NDAT];
    for e in data.iter_mut() {
        *e = rng.uni();
    }

    let mut bdata = [0.0f64; NBINS];
    bin_it(&data, NBINS, &mut bdata);
    bdata.sort_by(|a, b| a.partial_cmp(b).unwrap());

    const HBINS: usize = 11;
    let bin_width = (xmax - xmin) / HBINS as f64;
    let mut histogram = [0; HBINS];

    for &value in &bdata {
        let bin_index = ((value - xmin) / bin_width).floor() as usize;
        if bin_index < HBINS {
            histogram[bin_index] += 1;
        }
    }
    let sum: usize = histogram.iter().sum();

    let scale_factor = 1.0 / (sum as f64 * bin_width);
    let histogram: Vec<f64> = histogram
        .iter()
        .map(|&count| count as f64 * scale_factor)
        .collect();

    let dx = (xmax - xmin) / histogram.len() as f64;
    let mut ymax = 0.0;
    let mut v = Vec::new();
    histogram.iter().enumerate().for_each(|(bin, &y)| {
        let x = xmin + dx * bin as f64;
        if y > ymax {
            ymax = y;
        }
        v.push((x, y));
        v.push((x + dx, y));
    });

    let dx = (xmax - xmin) / (10 * histogram.len()) as f64;
    let v2: Vec<(f64, f64)> = (0..10 * histogram.len())
        .map(|i| xmin + dx * i as f64)
        .map(|x| (x, gaussian_pdf(x, 0.5, sdv)))
        .collect();

    let mut ymax = 0.0;
    for t in [&v, &v2] {
        for (_, y) in t {
            if *y > ymax {
                ymax = *y;
            }
        }
    }

    plot(
        "fig1_9.png",
        vec![
            (String::from("Binned data"), v),
            (String::from("Normal distribution"), v2),
        ],
        xmin,
        xmax,
        ymax,
    );
}
