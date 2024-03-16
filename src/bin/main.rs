use marsaglia_rs::Marsaglia;
use std::collections::HashMap;
use ndarray::Array1;
use statrs::function::erf::erf;


fn frequency_test() {
    let p = 0.4;
    println!("\nFrequency of samples in [0;{p})");
    println!("k, #samples, frequency, error");
    for k in 1..12 {
        let m = 2_u64.pow(2 * k - 1);
        let mut rng = Marsaglia::new(12, 34, 56, 78);
        let n = (0..m).filter(|_| rng.uni() < p).count();
        let r = (n as f64) / (m as f64);
        println!("{k:3} {m:8} {r:.20} {:+e}", (r - p).abs());
    }
}
fn uniform_histogram() {
    let mut ymax = 0.0;
    let mut graphs = Vec::new();
    for m in [100, 10000] {
        let mut rng = Marsaglia::new(12, 34, 56, 78);
        let mut histogram: HashMap<usize, usize> = HashMap::new();

        (0..m).map(|_| (10.0 * rng.uni()) as usize).for_each(|bin| {
            let count = histogram.entry(bin).or_insert(0);
            *count += 1
        });

        println!("\nHistogram - {m} samples, 10 bins");
        let mut bins: Vec<usize> = histogram.keys().cloned().collect();
        bins.sort();
        let width = (m as f64).ln() as usize;

        for bin in &bins {
            println!("Bin {bin}: {:width$}", histogram[&bin]);
        }

        let dx = 1.0 / bins.len() as f64;
        let sum: f64 = dx * histogram.values().sum::<usize>() as f64;

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
        graphs.push((format!("{m}"), v));
    }
    plot0("plot0.png", graphs, ymax + 0.1);
}

fn min_max_test() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    let mut f1 = rng.uni();
    let mut c1 = 1;
    let mut f2 = f1;
    let mut c2 = 1;
    const M: u64 = 10_000_000_000;
    println!("\nmin, max and their frequencies over {M} samples");

    for _ in 0..M - 1 {
        let f = rng.uni();
        if f == f1 {
            c1 += 1;
        } else if f < f1 {
            f1 = f;
            c1 = 0;
        }
        if f == f2 {
            c2 += 1;
        } else if f > f2 {
            f2 = f;
            c2 = 0;
        }
    }
    println!("min {f1}; count {c1}");
    println!("max {f2}; count {c2}");
}

fn gaussian_histogram() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    let m = 10_000;
    let mut histogram: HashMap<i32, usize> = HashMap::new();
    (0..m)
        .map(|_| rng.gauss())
        .map(|z| (10.0 * z).floor() as i32)
        .filter(|&b| b >= -30 && b <= 30)
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
    plot1(
        "plot1.png",
        vec![(String::from("Histogram"), v)],
        ymax + 0.1,
    );
}

fn main() {
    frequency_test();

    uniform_histogram();
    min_max_test();

    gaussian_histogram();

    plot2(); //cdf
}

use plotters::prelude::*;

// reproduces fig 1.2 p. 14
fn plot0(fname: &str, graphs: Vec<(String, Vec<(f64, f64)>)>, ymax: f64) {
    let root_area = BitMapBackend::new(fname, (600, 400)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption("Normalised Histograms", ("sans-serif", 40))
        .build_cartesian_2d(0f64..1f64, 0f64..ymax)
        .unwrap();

    ctx.configure_mesh().x_desc("x").y_desc("H").draw().unwrap();

    graphs.into_iter().enumerate().for_each(|(i, (label, h))| {
        let colour = match i {
            0 => RED,
            _ => BLUE,
        };
        ctx.draw_series(
            AreaSeries::new(
                h.into_iter(),
                0.0,              // Baseline
                &colour.mix(0.2), // Make the series opac
            )
            .border_style(&colour), // Make a brighter border
        )
        .unwrap()
        .label(label)
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &colour));
    });

    ctx.configure_series_labels()
        .border_style(&BLACK)
        .background_style(&WHITE.mix(0.8))
        .position(SeriesLabelPosition::UpperRight)
        .margin(20)
        .draw()
        .unwrap();
}

use std::f64::consts::PI;

fn gaussian_pdf(x: f64, mean: f64, std_dev: f64) -> f64 {
    let exponent = -((x - mean) * (x - mean)) / (2.0 * std_dev * std_dev);
    let coefficient = 1.0 / (std_dev * (2.0 * PI).sqrt());
    coefficient * exponent.exp()
}

fn cauchy_pdf(x: f64, x0: f64, gamma: f64) -> f64 {
    1.0 / (PI * gamma * (1.0 + ((x - x0) / gamma).powi(2)))
}

fn cauchy_cdf(x: f64, x0: f64, gamma: f64) -> f64 {
    1.0 / std::f64::consts::PI * ((x - x0) / gamma).atan() + 0.5
}

fn uniform_pdf(x: f64, x0: f64, x1: f64) -> f64 {
    if x >= x0 && x <= x1 {
        1.0 / (x1 - x0)
    } else {
        0.0
    }
}

fn uniform_cdf(x: f64, x0: f64, x1: f64) -> f64 {
    if x < x0 {
        0.0
    } else if x < x1 {
        (x - x0) / (x1 - x0)
    } else {
        1.0
    }
}

// reporduces fig1.4 p.25
fn plot1(fname: &str, graphs: Vec<(String, Vec<(f64, f64)>)>, ymax: f64) {
    let root_area = BitMapBackend::new(fname, (600, 400)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption("Probability Densities", ("sans-serif", 40))
        .build_cartesian_2d(-3f64..3f64, 0f64..ymax)
        .unwrap();

    ctx.configure_mesh().x_desc("x").y_desc("f").draw().unwrap();

    graphs.into_iter().enumerate().for_each(|(i, (label, h))| {
        let colour = match i {
            0 => RED,
            _ => BLUE,
        };
        ctx.draw_series(
            AreaSeries::new(
                h.into_iter(),
                0.0,              // Baseline
                &colour.mix(0.2), // Make the series opac
            )
            .border_style(&colour), // Make a brighter border
        )
        .unwrap()
        .label(label)
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &colour));
    });

    let data = (-30..=30)
        .map(|x| x as f64 / 10.0)
        .map(|x| (x, gaussian_pdf(x, 0.0, 1.0)));
    //let v: Vec<(f64,f64)> = data.clone().collect();
    //println!("{:?}", v);
    //let colour=GREEN;
    let colour = Palette99::pick(2).mix(0.9);
    ctx.draw_series(LineSeries::new(data, &colour))
        .unwrap()
        .label("Gaussian")
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &colour));

    let data = (-30..=30)
        .map(|x| x as f64 / 10.0)
        .map(|x| (x, cauchy_pdf(x, 0.0, 1.0)));
    let colour = Palette99::pick(3).mix(0.9);
    ctx.draw_series(LineSeries::new(data, &colour))
        .unwrap()
        .label("Cauchy")
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &colour));

    let data = (-30..=30)
        .map(|x| x as f64 / 10.0)
        .map(|x| (x, uniform_pdf(x, -1.0, 1.0)));
    let colour = Palette99::pick(4).mix(0.9);
    ctx.draw_series(LineSeries::new(data, &colour))
        .unwrap()
        .label("Uniform")
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &colour));

    ctx.configure_series_labels()
        .border_style(&BLACK)
        .background_style(&WHITE.mix(0.8))
        .position(SeriesLabelPosition::UpperRight)
        .margin(20)
        .draw()
        .unwrap();
}

// reproduces fig 1.5 p. 26
fn plot2() {
    // Generate x-values
    let x: Array1<f64> = Array1::linspace(-3.0, 3.0, 1000);

    // Calculate CDFs
    let gauss_cdf = x.map(|x| 0.5 * (1.0 + erf(x / (2.0f64.sqrt()))));
    let cau_cdf = x.map(|x| cauchy_cdf(*x, 0.0,1.0));
    let uni_cdf = x.map(|x| uniform_cdf(*x, -1.0,1.0));

    // Set up the plot
    let root = BitMapBackend::new("plot2.png", (800, 600)).into_drawing_area();
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
            &RED,
        ))
        .unwrap()
        .label("Gaussian CDF")
        .legend(|(x, y)| Path::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(
            x.iter().cloned().zip(cau_cdf.iter().cloned()),
            &BLUE,
        ))
        .unwrap()
        .label("Cauchy CDF")
        .legend(|(x, y)| Path::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .draw_series(LineSeries::new(
            x.iter().cloned().zip(uni_cdf.iter().cloned()),
            &GREEN,
        ))
        .unwrap()
        .label("Uniform CDF")
        .legend(|(x, y)| Path::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();
}
