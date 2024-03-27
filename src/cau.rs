pub fn cauchy_pdf(x: f64, x0: f64, gamma: f64) -> f64 {
    1.0 / (std::f64::consts::PI * gamma * (1.0 + ((x - x0) / gamma).powi(2)))
}

pub fn cauchy_cdf(x: f64, x0: f64, gamma: f64) -> f64 {
    1.0 / std::f64::consts::PI * ((x - x0) / gamma).atan() + 0.5
}

// CAUCHY CUMULATIVE DISTRIBUTIO FUNCTION.
pub fn cau_df(x: f64) -> f64 {
    0.5 + x.atan() / std::f64::consts::PI
}

pub fn cau_qdf(x: f64) -> f64 {
    let y = cau_df(x);
    if y > 0.5 {
        1.0 - y
    } else {
        y
    }
}

pub fn cau_xq(q: f64) -> f64 {
    let qq = std::f64::consts::PI * (q - 0.5);
    qq.tan()
}
