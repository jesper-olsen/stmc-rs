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
    let exponent = -((x - mean) * (x - mean)) / (2.0 * std_dev * std_dev);
    let coefficient = 1.0 / (std_dev * (2.0 * PI).sqrt());
    coefficient * exponent.exp()
}

/// GAUSSIAN CUMULATIVE DISTRIBUTION FUNCTION.
pub fn gau_df(x: f64) -> f64 {
    //0.5 *(1.0 + erf(x / (2.0f64.sqrt())))
    0.5 * (1.0 + gamma::error_f(x / (2.0f64.sqrt())))
}
