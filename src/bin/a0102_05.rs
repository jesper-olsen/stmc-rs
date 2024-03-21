use marsaglia_rs::marsaglia::Marsaglia;

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
            c1 = 1;
        }
        if f == f2 {
            c2 += 1;
        } else if f > f2 {
            f2 = f;
            c2 = 1;
        }
    }
    println!("min {f1}; count {c1}");
    println!("max {f2}; count {c2}");
}

fn main() {
    min_max_test();
}
