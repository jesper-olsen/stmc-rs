use stmc_rs::marsaglia::Marsaglia;
use stmc_rs::steb::steb1;

fn main() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    const KMAX: usize = 10;
    let mut data = [0.0f64; KMAX];
    let mut w = [0.0f64; KMAX];
    for k in 1..=KMAX {
        let ndat = if k == KMAX { 2 } else { 2usize.pow(k as u32) };
        data[k - 1] = (0..ndat).map(|_| rng.uni()).sum::<f64>() / ndat as f64;
        w[k - 1] = ndat as f64;
        println!(" K,NDAT,DATA(K): {k:3}, {ndat:3}, {:2.4}", data[k - 1]);
    }
    let (xm, xe) = steb1(&data, &mut w);
    println!("UNIFORM RANDOM NUMBERS: XM={xm} +/- {xe}");
}
