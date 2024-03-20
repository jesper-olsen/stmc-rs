use marsaglia_rs::Marsaglia;

fn steb0(data: &[f64]) -> (f64, f64) {
    let xm = data.iter().sum::<f64>() / data.len() as f64;
    let xv = data.iter().map(|x| (x - xm).powi(2)).sum::<f64>() / (data.len() as f64 - 1.0);
    let xe = (xv / (data.len() as f64)).sqrt();
    (xm, xe)
}

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
    const KMAX: usize = 6;
    const NDT0: usize = 2usize.pow(KMAX as u32);
    const NRPT: usize = 10_000;
    let mut dat0 = [0.0; NDT0];
    let mut wght = [0.0; KMAX];
    let mut va0 = [0.0; NRPT];
    let mut va1 = [0.0; NRPT];
    let mut eb0 = [0.0; NRPT];
    let mut eb1 = [0.0; NRPT];
    println!("\nKMAX,NRPT: {KMAX},{NRPT}");
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    for j in 0..NRPT {
        let mut i = 0;
        let mut dat1 = [0.0; KMAX];
        for k in 0..KMAX {
            let ndat = if k == KMAX - 1 {
                2
            } else {
                2usize.pow(k as u32 + 1)
            };
            for _ in 0..ndat {
                dat0[i] = rng.uni();
                dat1[k] += dat0[i];
                i += 1;
            }
            wght[k] = ndat as f64;
            dat1[k] /= wght[k];
        }
        if i != NDT0 {
            panic!("ERROR IN COUNTING i=1,...,NDT0.");
        }
        let (xm1, xe1) = steb1(&dat1, &mut wght);
        let (xm0, xe0) = steb0(&dat0);
        va1[j] = xe1.powi(2);
        va0[j] = xe0.powi(2);
        eb1[j] = xe1;
        eb0[j] = xe0;
        if j == NRPT - 1 {
            println!("UNIFORM RANDOM NUMBERS: XM1 = {xm1:.6} +/- {xe1:.6}");
            println!("UNIFORM RANDOM NUMBERS: XM0 = {xm0:.6} +/- {xe0:0.6}");
        }
    }
    let (eb0m, xe0) = steb0(&eb0);
    let (eb1m, xe1) = steb0(&eb1);
    let (xv0m, ve0) = steb0(&va0);
    let (xv1m, ve1) = steb0(&va1);
    println!("\nVARIANCE                     ERROR BAR");
    println!(" steb1: {xv1m:.6} +/- {ve1:.6}, {eb1m:.6} +/- {xe1:.6}");
    println!(" steb0: {xv0m:.6} +/- {ve0:.6}, {eb0m:.6} +/- {xe0:.6}");
}
