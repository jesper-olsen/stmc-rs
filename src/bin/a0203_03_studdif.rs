use stmc_rs::student::stud_df;

// Table 2.8 p.80
const DATA: [(usize, usize, f64); 9] = [
    (512, 512, 0.048),
    (64, 64, 0.0054),
    (16, 16, 0.0083),
    (16, 8, 0.020),
    (8, 8, 0.023),
    (16, 4, 0.072),
    (4, 4, 0.030),
    (3, 3, 0.047),
    (2, 2, 0.11),
];

const XM1: f64 = 1.0;
const XM2: f64 = 1.2;
const EB1: f64 = 0.05f64;
const EB2: f64 = 0.05f64;

///STUDENT'S DIFFERENCE TEST: COMPARISION OF TWO MEANS.
fn main() {
    for (n, m, q) in DATA.iter() {
        let nf = n + m - 2; // degrees of freedom
        let n = *n as f64;
        let m = *m as f64;
        let s12 = n * EB1.powi(2); // squared variance
        let s22 = n * EB2.powi(2);
        let s2 = ((n - 1.0) * s12 + (m - 1.0) * s22) / nf as f64;
        let s2 = (1.0 / n + 1.0 / m) * s2;
        let t = (XM1 - XM2).abs() / s2.sqrt();
        let qa = 2.0 * stud_df(-t, nf);
        println!("N,M,S2: {n:3}, {m:3}, {s2:1.4} Q: {q:1.4} vs {qa}",);
    }
}
