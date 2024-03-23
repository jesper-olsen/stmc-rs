use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::qtiles;
use marsaglia_rs::steb::steb0;

fn fractiles(q: f64) -> (f64, f64) {
    const NDAT: usize = 10000;
    const ND2: usize = 2;
    let mut var = [0.0f64; NDAT];
    let mut data = [0.0f64; NDAT];
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    for e in var.iter_mut() {
        data[0] = rng.gauss();
        data[1] = rng.gauss();
        let (_, v) = steb0(&data[0..ND2]);
        *e = 2.0 * v.powi(2);
    }
    var.sort_by(|a, b| a.partial_cmp(b).unwrap());
    qtiles(&var, q).unwrap()
}

/// Variance fractiles: experimental verification.
fn main() {
    let q = 0.15;
    let (vq1, vq2) = fractiles(q);
    println!("q: {q} vq1: {vq1} vq2: {vq2}");
}

#[cfg(test)]
mod tests {
    use crate::fractiles;
    #[test]
    fn test_fractiles() {
        let (vq1, vq2) = fractiles(0.15);
        assert_eq!(vq1, 0.03470057507208726);
        assert_eq!(vq2, 2.045777357830865);
    }
}
