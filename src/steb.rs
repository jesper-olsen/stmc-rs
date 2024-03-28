///C  INPUT:   ARRAY X(N) WITH (GAUSSIAN) DATA.
///C  OUTPUT:  MEAN VALUE XM, UNBIASED VARIANCE XV AND ERROR BAR XE.
pub fn steb0(data: &[f64]) -> (f64, f64) {
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

///C  INPUT:   ARRAY X (GAUSSIAN) DATA AND WEIGHT FACTORS.
///C  OUTPUT:  MEAN VALUE XM, AND ITS ERROR BAR XE.
///C           WEIGHT FACTORS NORMALIZED TO ONE.
pub fn steb1(data: &[f64], w: &mut [f64]) -> (f64, f64) {
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

///C  INPUT:   DATA (GAUSSIAN), ERROR BARS  AND  (OPTIONAL) WEIGHT
///C           FACTORS. (IF WEIGHT FACTORS ARE NOT GIVEN: PUT
///C           W(1) < 0  AND WEIGHT FACTORS WILL WE CALCULATED
///C           FROM THE ERROR BARS EB.)
///C  OUTPUT:  MEAN VALUE XM AND ITS ERROR BAR XE.
///C           WEIGHT FACTORS ARE RETURNED NORMALIZED TO ONE.
pub fn steb2(data: &[f64], eb: &[f64], w: &mut [f64]) -> (f64, f64) {
    if w[0] < 0.0 {
        for (w, e) in w.iter_mut().zip(eb.iter()) {
            *w = 1.0 / (*e).powi(2);
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

///CALCULATION OF JACKKNIFE ESTIMATORS FROM  N  JACKKNIFE BINS  FJ(N):
///- THE MEAN VALUE:   FM,
///- THE VARIANCE:     FV,
///- THE ERROR BAR:    FE=SQRT(FV/N), (STANDARD DEVIATION OF FM).   
pub fn stebj0(fj: &[f64]) -> (f64, f64, f64) {
    assert!(fj.len() > 1, "stebj0: fj too short");
    let fm = fj.iter().sum::<f64>() / fj.len() as f64;
    let fv = fj.iter().map(|x| (x - fm).powi(2)).sum::<f64>() * (fj.len() - 1) as f64;
    let fe = (fv / fj.len() as f64).sqrt();
    (fm, fv, fe)
}

///   PURPOSE: CALCULATION OF THE JACKKNIFE ESTIMATORS FJM  AND
///            THE UNBIASED ESTIMATOR FMM FOR  THE MEAN VALUE <F>,
///            JACKKNIFE ESTIMATOR FV FOR THE VARIANCE <F*F>-<F><F>,
///            AND ERROR BAR, DEFINED AS STANDARD DEVIATION OF FJM
///            FE  := SQRT(FMV) = SQRT(FV/N).   
///                                                                    
///   INPUT:
///                                                                    
///   FJ      VECTOR CONTAINING THE MEASUREMENTS TO BE AVERAGED;    
///   FMEAN   MEAN VALUE FROM COMPLETE STATISTICS.
///
///   OUTPUT:
///
///   FJM     JACKKNIFE  MEAN VALUE;
///   FMM     UNBIASED MEAN VALUES;
///   FV      VARIANCE;            
///   FE      ERROR BAR.
pub fn stebj1(fj: &[f64], fmean: f64) -> (f64, f64, f64, f64) {
    assert!(fj.len() >= 2, "stebj1: fj too small!");
    let (fjm, fv, fe) = stebj0(fj);
    let (fmm, biasm) = bias(fj, fmean);
    (fjm, fmm, fv, fe)
}

//  ESTIMATES THE BIAS (ASSUMPTION: BIAS=CONST/N)  BIASM  AND
// THE BIAS-CORRECTED MEAN VALUE  XMM  FROM THE FOLLOWING
// INPUT:  XJ(N) JACKKNIFE BINS AND
//         XM    MEAN VALUE FROM ENTIRE STATISTICS.
fn bias(xj: &[f64], xm: f64) -> (f64, f64) {
    let xjm = xj.iter().sum::<f64>() / xj.len() as f64;
    let biasm = (xj.len() - 1) as f64 * (xm - xjm);
    (xm + biasm, biasm)
}
