use marsaglia_rs::marsaglia::Marsaglia;

fn steb1(data: &[f64], w: &mut [f64]) -> (f64, f64) {
    //C  INPUT:   ARRAY X (GAUSSIAN) DATA AND WEIGHT FACTORS.
    //C  OUTPUT:  MEAN VALUE XM, AND ITS ERROR BAR XE.
    //C           WEIGHT FACTORS NORMALIZED TO ONE.
    assert_eq!(data.len(), w.len());
    let wnorm = w.iter().sum::<f64>();
    for e in w.iter_mut() {
        *e /= wnorm;
    }
    let xm = data
        .iter()
        .zip(w.iter())
        .fold(0.0, |acc, (&x, &y)| acc + x * y);
    let xv = data
        .iter()
        .zip(w.iter())
        .fold(0.0, |acc, (&x, &y)| acc + y * (x - xm).powi(2));
    let xe = (xv / (data.len() - 1) as f64).sqrt();
    (xm, xe)
}

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
