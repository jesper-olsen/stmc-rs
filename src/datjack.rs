// CALCULATION OF  N  JACKKNIFE BINS  XJ(N)  FOR  N  DATA IN  X(N).
pub fn datjack(x: &[f64], xj: &mut [f64]) {
    let xsum = x.iter().sum::<f64>();
    for (e, f) in xj.iter_mut().zip(x.iter()) {
        *e = (xsum - *f) / (x.len() - 1) as f64;
    }
    //for i in 0..x.len() {
    //    xj[i]=(xsum-x[i])/(x.len()-1) as f64;
    //}
}
