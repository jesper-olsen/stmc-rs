use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::plot::plot;
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

    let xmin = 0.0;
    let xmax = 1.0;
    plot::plot("fig1_2.png", "Normalised Histograms", "x", "H", graphs, xmin, xmax, ymax);
}

fn main() {
    uniform_histogram();
}
