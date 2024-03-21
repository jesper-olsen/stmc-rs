use std::io;

pub mod marsaglia;
pub mod gamma;

fn get_input() -> String {
    let mut s = String::new();
    io::stdin().read_line(&mut s).expect("Failed to read line");
    String::from(s.trim())
}

pub fn get_float(msg: &str) -> f64 {
    loop {
        println!("{msg}");
        let input: String = get_input();
        if let Ok(num) = input.trim().parse::<f64>() {
            return num;
        } else {
            println!("Invalid input. Please enter a valid floating-point number.");
        };
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
