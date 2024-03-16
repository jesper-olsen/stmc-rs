use marsaglia_rs::Marsaglia;
use std::collections::HashMap;

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

        let sum: f64 = histogram.values().sum::<usize>() as f64;

        let dx = 1.0 / bins.len() as f64;
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
    plot0("histogram0.png", graphs, ymax + 0.1);
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
    let m = 1000000;
    let mut histogram: HashMap<i32, usize> = HashMap::new();
    (0..m)
        .map(|_| rng.gauss())
        .map(|z| (10.0 * z) as i32)
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

//    let sum: f64 = histogram.values().sum::<usize>() as f64;
//    let mut mx = 0.0;
//    let dx = 1.0 / bins.len() as f64;
//    let mut v = Vec::new();
//    for bin in bins {
//        let x = bin as f64 / 10.0;
//        let y = histogram[&bin] as f64 / sum;
//        if y > mx {
//            mx = y;
//        }
//        v.push((x, y));
//        v.push((x + dx, y));
//    }
//    v
}

fn main() {
    frequency_test();

    uniform_histogram();
    min_max_test();

    gaussian_histogram();

}

use plotters::prelude::*;

fn plot0(fname: &str, graphs: Vec<(String, Vec<(f64, f64)>)>, ymax: f64) {
    let root_area = BitMapBackend::new(fname, (600, 400)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption("Histogram", ("sans-serif", 40))
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
