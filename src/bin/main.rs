

use marsaglia_rs::Marsaglia;

fn main() {
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

    use std::collections::HashMap;
    for m in [100, 100000] {
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
        for bin in bins {
            println!("Bin {bin}: {:width$}", histogram[&bin]);
        }
    }

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
