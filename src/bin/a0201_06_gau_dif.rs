use marsaglia_rs::gamma::{error_f, gamma_ln, gamma_p};
use marsaglia_rs::get_float;
use std::io;

fn gaudif(xm1: f64, eb1: f64, xm2: f64, eb2: f64) -> f64 {
    let sigma = (eb1.powi(2) + eb2.powi(2)).sqrt();
    let xx = (xm1 - xm2).abs() / (sigma * 2.0f64.sqrt());
    1.0 - error_f(xx)
}

fn main() {
    println!("DIFFERENCE TEST: COMPARISON OF TWO MEANS.");
    println!("INPUT: TWO MEAN VALUES AND THEIR ERROR BARS (NOT VARIANCES!");
    println!("LIKELIHOOD FOR THE OBSERVED DISCREPANCE TO BE DUE TO CHANCE");
    let xm1 = get_float("ENTER 1. MEAN VALUE (DATA POINT):");
    let eb1 = get_float("ENTER ERROR BAR OF 1. DATA POINT:");
    let xm2 = get_float("ENTER 2. MEAN VALUE (DATA POINT):");
    let eb2 = get_float("ENTER ERROR BAR OF 2. DATA POINT:");

    let q = gaudif(xm1, eb1, xm2, eb2);

    println!("Q = {q}");
    println!("OR EQUIVALENTLY: {:.2}%", q * 100.0);
}
