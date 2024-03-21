use marsaglia_rs::gamma::{error_f};

// GAUSSIAN: CONFIDENCE LIMIT FOR MULTIPLE STANDARD DEVIATIONS.

fn main() {
    const N: usize = 6;
    let mut p = [0.0; N];
    let mut q = [0.0; N];
    let mut na = [0; N];

    for i in 0..N {
        na[i] = i + 1;
        let s = (i + 1) as f64;
        let x = s / f64::sqrt(2.0);
        p[i] = error_f(x);
        q[i] = 0.5 * (1.0 - p[i]);
    }

    println!(
        "     n{:>11}{:>11}{:>11}{:>11}{:>11}{:>11}",
        na[0], na[1], na[2], na[3], na[4], na[5]
    );
    println!(
        "     p{:11.2}{:11.2}{:11.2}{:11.2}{:11.2}{:11.2}",
        p[0], p[1], p[2], p[3], p[4], p[5]
    );
    println!(
        "     q{:11.2e}{:11.2e}{:11.2e}{:11.2e}{:11.2e}{:11.2e}",
        q[0], q[1], q[2], q[3], q[4], q[5]
    );
}
