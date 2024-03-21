use marsaglia_rs::marsaglia::Marsaglia;
use std::f64::consts::PI;

fn main() {
    println!("MOMENTS 1-4 FOR UNIFORM RANDOM NUMBERS AND 1-2 FOR CAUCHY RANDOM NUMERS:\n");
    println!(" K    ndat       U1       U2       U3        U4        C1         C2\n");

    for k in 1..=10 {
        let ndat: usize = 1 << (2 * k - 1);
        let mut xmcau = 0.0;
        let mut x2cau = 0.0;
        let mut xmo = [0.0; 4];
        let mut xmom = [0.0; 4];

        let mut rng = Marsaglia::new(12, 34, 56, 78);
        for _ in 0..ndat {
            let xr: f64 = rng.uni();
            xmcau += (2.0 * PI * xr).tan();
            x2cau += (2.0 * PI * xr).tan().powi(2);
            xmo[0] = xr;
            xmom[0] += xr;
            for i in 1..4 {
                xmo[i] = xmo[i - 1] * xr;
                xmom[i] += xmo[i];
            }
        }

        // Normalization
        xmcau /= ndat as f64;
        x2cau /= ndat as f64;
        for e in &mut xmom {
            *e /= ndat as f64;
        }

        println!(
            "{k:>2} {ndat:>8} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {xmcau:>10.4} {x2cau:>12.1}",
            xmom[0], xmom[1], xmom[2], xmom[3]
        );
    }
}
