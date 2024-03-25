use marsaglia_rs::f::f_xq;

// Confidence intervals for F ratio estimates
// Reproduces Table 2.12 p.88
fn calc_table() -> Vec<f64> {
    let kdat = 10;
    let pc1 = 0.7;
    let pc2 = 0.95;

    let q1 = 0.5 * (1.0 - pc1);
    let p1 = 1.0 - q1;
    let q2 = 0.5 * (1.0 - pc2);
    let p2 = 1.0 - q2;
    let nda1 = 2u32.pow(kdat);
    let nf1 = nda1 - 1;

    println!("\n NDA1 = {nda1}:\n");
    println!("   SOME CONFIDENCE INTERVALS FOR F RATIO ESTIMATES.\n");
    println!("                             q               q             1-q             1-q\n");
    println!("          NDA2            .025            .150            .850            .975");
    let mut v = Vec::new();
    for k2 in -4..=4isize {
        let expo = kdat as isize + k2;
        let nda2 = 2u32.pow(expo as u32) as u32;
        let nf2 = nda2 - 1;
        let fdo2 = 1.0 / f_xq(p2, nf1, nf2);
        let fdo1 = 1.0 / f_xq(p1, nf1, nf2);
        let fup1 = 1.0 / f_xq(q1, nf1, nf2);
        let fup2 = 1.0 / f_xq(q2, nf1, nf2);
        println!(
            "   {nda2:11}      {fdo2:10.3}      {fdo1:10.3}      {fup1:10.3}      {fup2:10.3}"
        );
        v.push(fdo2);
        v.push(fdo1);
        v.push(fup1);
        v.push(fup2);
    }
    v
}

fn main() {
    calc_table();
}

#[cfg(test)]
mod tests {
    use crate::calc_table;

    fn round_to_three_decimals(value: f64) -> f64 {
        (value * 1000.0).round() / 1000.0
    }

    #[test]
    fn f_xq_test() {
        let r: Vec<f64> = vec![
            0.676, 0.813, 1.192, 1.396, 0.759, 0.865, 1.140, 1.282, 0.819, 0.900, 1.105, 1.208,
            0.859, 0.923, 1.082, 1.160, 0.885, 0.937, 1.067, 1.130, 0.900, 0.946, 1.058, 1.113,
            0.909, 0.951, 1.053, 1.103, 0.914, 0.953, 1.051, 1.098, 0.916, 0.955, 1.049, 1.095,
        ];

        let r2: Vec<f64> = calc_table()
            .iter()
            .map(|x| round_to_three_decimals(*x))
            .collect();
        assert_eq!(r, r2);
    }
}
