use marsaglia_rs::Marsaglia;
use plotters::prelude::*;
use std::collections::HashMap;

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
    plot("fig1_2.png", graphs, ymax + 0.1);
}

fn main() {
    uniform_histogram();
}

// reproduces fig 1.2 p. 14
fn plot(fname: &str, graphs: Vec<(String, Vec<(f64, f64)>)>, ymax: f64) {
    println!("Saving plot: {fname}"); 
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
