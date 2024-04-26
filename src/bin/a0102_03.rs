use stmc_rs::marsaglia::Marsaglia;

fn frequency_test() {
    let p = 0.4;
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    println!("\nSampling; Frequency of samples in [0;{p})");
    for k in 1..12 {
        let m = 2_u64.pow(2 * k - 1);
        let n = (0..m).filter(|_| rng.uni() < p).count();
        let r = (n as f64) / (m as f64);
        println!("K,N,R,ERROR: {k:3} {m:8} {r:.20} {:+e}", (r - p).abs());
    }
}

fn main() {
    frequency_test();
}
