use marsaglia_rs::Marsaglia;

fn steb0(data: &[f64]) -> (f64, f64) {
    let xm = data.iter().sum::<f64>() / data.len() as f64;
    let xv = data.iter().map(|x| (x - xm).powi(2)).sum::<f64>() / (data.len() as f64 - 1.0);
    let xe = (xv / (data.len() as f64)).sqrt();
    (xm, xe)

    // 1-pass
    //    let mut xm = 0.0f64;
    //    let mut xv = 0.0f64;
    //    let mut n: usize = 1;
    //    for &x in data {
    //        n += 1;
    //        let delta = x - xm;
    //        xm += delta / n as f64;
    //        xv += delta * (x - xm);
    //    }
    //    if n > 1 {
    //        xv /= (n - 1) as f64; // Unbiased estimator for sample variance
    //    } else {
    //        xv = f64::NAN; // Handle single-element input
    //    }
    //
    //    let xe = (xv / (data.len() as f64)).sqrt();
    //    return (xm, xe);
}

fn main() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);

    let mut data = [0.0f64; 1024];
    for e in data.iter_mut() {
        *e = rng.uni();
    }
    let (xm, xe) = steb0(&data);
    println!("UNIFORM  RANDOM NUMBERS: {xm} {xe}");

    let mut data = [0.0f64; 2048];
    for e in data.iter_mut() {
        *e = rng.gauss();
    }
    let (xm, xe) = steb0(&data);
    println!("GAUSSIAN RANDOM NUMBERS: {xm} {xe}");

    let mut data = [0.0f64; 1024];
    for e in data.iter_mut() {
        *e = rng.cauchy();
    }
    let (xm, xe) = steb0(&data);
    println!("CAUCHY   RANDOM NUMBERS: {xm} {xe}");
}
