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

fn steb2(data: &[f64], eb: &[f64], w: &mut [f64]) -> (f64, f64) {
    //C  INPUT:   DATA (GAUSSIAN), ERROR BARS  AND  (OPTIONAL) WEIGHT
    //C           FACTORS. (IF WEIGHT FACTORS ARE NOT GIVEN: PUT
    //C           W(1) < 0  AND WEIGHT FACTORS WILL WE CALCULATED
    //C           FROM THE ERROR BARS EB.)
    //C  OUTPUT:  MEAN VALUE XM AND ITS ERROR BAR XE.
    //C           WEIGHT FACTORS ARE RETURNED NORMALIZED TO ONE.

    if w[0] < 0.0 {
        for i in 0..w.len() {
            w[i] = 1.0 / eb[i].powi(2);
        }
    }

    let wnorm: f64 = w.iter().sum();
    for e in w.iter_mut() {
        *e /= wnorm;
    }

    //weighted mean
    let xm = data
        .iter()
        .zip(w.iter())
        .fold(0.0, |acc, (&x, &y)| acc + x * y);

    //weighted error
    let xe = eb
        .iter()
        .zip(w.iter())
        .fold(0.0, |acc, (&x, &y)| acc + x.powi(2) * y.powi(2))
        .sqrt();

    //variance
    //xe*f64::sqrt(data.len() as f64
    (xm, xe)
}

fn main() {
    const KMAX: usize = 10;
    const NMAX: usize = 2usize.pow(KMAX as u32 - 1);
    let mut data = [0.0f64; KMAX];
    let mut dat1 = [0.0f64; NMAX];
    let mut eb = [0.0f64; KMAX];
    let mut wght = [0.0f64; KMAX];
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    wght[0] = -1.0;
    for k in 0..KMAX {
        let ndat = if k == KMAX - 1 {
            2
        } else {
            2usize.pow(1 + k as u32)
        };
        let a = &mut dat1[0..ndat];
        for e in a.iter_mut() {
            *e = rng.uni();
        }
        (data[k], eb[k]) = steb0(a);
    }
    let (xm, xe) = steb2(&data, &eb, &mut wght);
    println!("UNIFORM RANDOM NUMBERS: XM={xm:.7} +/- {xe:.7}");
}
