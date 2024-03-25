use std::f64::consts::PI;
use std::io::{self, Write};

pub mod beta;
pub mod chi2;
pub mod gamma;
pub mod gau;
pub mod marsaglia;
pub mod plot;
pub mod steb;
pub mod student;
pub mod f;

fn get_input() -> String {
    let mut s = String::new();
    io::stdin().read_line(&mut s).expect("Failed to read line");
    String::from(s.trim())
}

pub fn get_number<T: std::str::FromStr>(msg: &str) -> T {
    loop {
        println!("{msg}");
        let input: String = get_input();
        if let Ok(num) = input.trim().parse::<T>() {
            return num;
        } else {
            println!("Invalid input. Please enter a valid number.");
        };
    }
}

pub fn yes(q: &str, y: &str, n: &str) -> bool {
    loop {
        print!("{q}\n** ");
        let _ = io::stdout().flush();
        if let Some(c) = get_input().chars().next() {
            match c {
                'Y' | 'y' => {
                    println!("{y}");
                    return true;
                }
                'N' | 'n' => {
                    println!("{n}");
                    return false;
                }
                _ => println!(" Please answer Yes or No."),
            }
        }
    }
}

pub fn qtiles(x: &[f64], q: f64) -> Option<(f64, f64)> {
    let n = x.len();
    let nq = (q * n as f64) as usize;

    if nq == 0 || nq >= n {
        return None; // Return None if NQ is invalid
    }

    let w2 = q * (n + 1) as f64 - nq as f64;
    let w1 = 1.0 - w2;

    let xq1 = w1 * x[nq - 1] + w2 * x[nq];
    let xq2 = w1 * x[n - nq] + w2 * x[n - nq - 1];

    Some((xq1, xq2))
}

pub fn cauchy_pdf(x: f64, x0: f64, gamma: f64) -> f64 {
    1.0 / (PI * gamma * (1.0 + ((x - x0) / gamma).powi(2)))
}

pub fn cauchy_cdf(x: f64, x0: f64, gamma: f64) -> f64 {
    1.0 / std::f64::consts::PI * ((x - x0) / gamma).atan() + 0.5
}

pub fn uniform_pdf(x: f64, x0: f64, x1: f64) -> f64 {
    if x >= x0 && x <= x1 {
        1.0 / (x1 - x0)
    } else {
        0.0
    }
}

pub fn uniform_cdf(x: f64, x0: f64, x1: f64) -> f64 {
    if x < x0 {
        0.0
    } else if x < x1 {
        (x - x0) / (x1 - x0)
    } else {
        1.0
    }
}
