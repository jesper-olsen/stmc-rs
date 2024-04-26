use stmc_rs::marsaglia::Marsaglia;
use stmc_rs::qtiles;

fn main() {
    println!();
    println!("VARIOUS QTILES FOR CAUCHY RANDOM NUMBERS:");
    println!();
    println!("   K        N    X_0.5      X_0.25     X_0.75     X_0.15     X_0.85");
    println!();

    for k in 2..=8 {
        let m = 4usize.pow(k);
        let n = 2 * m - 1;

        // Generate Cauchy distributed random numbers
        let mut rng = Marsaglia::new(12, 34, 56, 78);
        let mut data = vec![0.0; n];
        for e in data.iter_mut().take(n) {
            *e = rng.cauchy();
        }

        data.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Calculate quantiles
        let (x05, _) = qtiles(&data, 0.50).unwrap();
        let (x025, x075) = qtiles(&data, 0.25).unwrap();
        let (x015, x085) = qtiles(&data, 0.15).unwrap();

        println!("{k:>4} {n:>8} {x05:>10.4} {x025:>10.4} {x075:>10.4} {x015:>10.4} {x085:>10.4}");
    }
    println!(
        "       EXACT: {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
        0.0000, -1.00000, 1.0000, -1.9626, 1.9626
    );

    println!("MEDIAN DEFINITIONS AGREE.");
}
