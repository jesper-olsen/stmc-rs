use marsaglia_rs::cau::{cauchy_cdf, cauchy_pdf};
use marsaglia_rs::gau::{gau_df, gau_pd};
use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::plot;
use marsaglia_rs::{uniform_cdf, uniform_pdf};
use ndarray::Array1;
use plotters::prelude::*;
use std::collections::HashMap;

fn gaussian_histogram() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    let m = 10_000;
    let mut histogram: HashMap<i32, usize> = HashMap::new();
    (0..m)
        .map(|_| rng.gauss())
        .map(|z| (10.0 * z).floor() as i32)
        .filter(|&b| (-30..=30).contains(&b))
        //.filter(|&b| b >= -30 && b <= 30)
        .for_each(|bin| {
            let count = histogram.entry(bin).or_insert(0);
            *count += 1
        });

    println!("\nHistogram - {} samples, 60 bins", 2 * m);
    let mut bins: Vec<i32> = histogram.keys().cloned().collect();
    bins.sort();
    let width = (2.0 * m as f64).ln() as usize;
    for bin in &bins {
        println!("Bin {bin:2}: {:width$}", histogram[&bin]);
    }

    let dx = 6.0 / bins.len() as f64;
    let sum: f64 = dx * histogram.values().sum::<usize>() as f64;

    let mut ymax = 0.0;
    let mut v = Vec::new();
    for bin in bins {
        let x = bin as f64 / 10.0;
        let y = histogram[&bin] as f64 / sum;
        if y > ymax {
            ymax = y;
        }
        v.push((x, y));
        v.push((x + dx, y));
    }

    let dgaussian: Vec<(f64, f64)> = (-30..=30)
        .map(|x| x as f64 / 10.0)
        .map(|x| (x, gau_pd(x, 0.0, 1.0)))
        .collect();

    let dcauchy: Vec<(f64, f64)> = (-30..=30)
        .map(|x| x as f64 / 10.0)
        .map(|x| (x, cauchy_pdf(x, 0.0, 1.0)))
        .collect();

    let duniform: Vec<(f64, f64)> = (-30..=30)
        .map(|x| x as f64 / 10.0)
        .map(|x| (x, uniform_pdf(x, -1.0, 1.0)))
        .collect();

    plot::plot(
        "fig1_4.png",
        "Histogram",
        "x",
        "f",
        vec![
            (String::from("Histogram CDF"), v.clone()),
            (String::from("Gaussian CDF"), dgaussian),
            (String::from("Cauchy CDF"), dcauchy),
            (String::from("Uniform"), duniform),
        ],
        -3.0..3.0,
        0.0..ymax,
    );
}

// reproduces fig 1.5 p. 26
fn plot2(fname: &str) {
    println!("Saving {fname}");
    // Generate x-values
    let x: Array1<f64> = Array1::linspace(-3.0, 3.0, 1000);

    // Calculate CDFs
    let gauss_cdf = x.map(|x| gau_df(*x));
    let cau_cdf = x.map(|x| cauchy_cdf(*x, 0.0, 1.0));
    let uni_cdf = x.map(|x| uniform_cdf(*x, -1.0, 1.0));

    // Set up the plot
    let root = BitMapBackend::new(fname, (800, 600)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "CDF Comparison: Gaussian vs Cauchy vs Uniform",
            ("sans-serif", 20).into_font(),
        )
        .margin(5)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-3.0..3.0, 0.0..1.0)
        .unwrap();

    // Plot CDFs
    chart
        .configure_mesh()
        .x_desc("X")
        .y_desc("F")
        .x_labels(5)
        .y_labels(10)
        .draw()
        .unwrap();

    chart
        .draw_series(LineSeries::new(
            x.iter().cloned().zip(gauss_cdf.iter().cloned()),
            RED,
        ))
        .unwrap()
        .label("Gaussian CDF")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));

    chart
        .draw_series(LineSeries::new(
            x.iter().cloned().zip(cau_cdf.iter().cloned()),
            BLUE,
        ))
        .unwrap()
        .label("Cauchy CDF")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));

    chart
        .draw_series(LineSeries::new(
            x.iter().cloned().zip(uni_cdf.iter().cloned()),
            GREEN,
        ))
        .unwrap()
        .label("Uniform CDF")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()
        .unwrap();
}

fn main() {
    gaussian_histogram();

    plot2("fig1_5.png");
}
