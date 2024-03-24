use std::f64::consts::PI;
use crate::gamma;

/// COMPARISION OF TWO MEANS: (GAUSSIAN DIFFERENCE TEST).
/// INPUT: TWO GAUSSIAN DATA POINTS, MEAN VALUES AND ERROR BARS.
/// OUPUT: LIKELIHOOD Q  THAT THE DISCREPANCY IS DUE TO CHANCE.
pub fn gaudif(xm1: f64, eb1: f64, xm2: f64, eb2: f64) -> f64 {
    let sigma = (eb1.powi(2) + eb2.powi(2)).sqrt();
    let xx = (xm1 - xm2).abs() / (sigma * 2.0f64.sqrt());
    1.0 - gamma::error_f(xx)
}

/// GAUSSIAN PROBABILITY DENSITY FUNCTION.
pub fn gau_pd(x: f64, mean: f64, std_dev: f64) -> f64 {
    (-0.5*(x-mean).powi(2)).exp()/(std_dev*(2.0*PI).sqrt())
}

/// GAUSSIAN CUMULATIVE DISTRIBUTION FUNCTION.
pub fn gau_df(x: f64) -> f64 {
    //0.5 *(1.0 + erf(x / (2.0f64.sqrt())))
    0.5 * (1.0 + gamma::error_f(x / (2.0f64.sqrt())))
}

// GAUSSIAN, PEAKED DISTRIBUTION FUNCTION. 
pub fn gau_qdf(x: f64) -> f64 {
    let f=0.5+0.5*gamma::error_f(x/2.0f64.sqrt());
    if f<=0.5 {
        f
    } else {
        1.0 - f
    } 
}
