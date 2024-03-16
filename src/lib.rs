// Rust implementation of the pseudo random number generator presented in
// TOWARD A UNIVERSAL RANDOM NUMBER GENERATOR
// George MARSAGLIA, Arif ZAMAN and Wai Wan TSANG
// Statistics & Probability Letters 8 (1990) 35-39

const N: usize = 97;
use std::f64::consts::PI;

pub struct Marsaglia {
    u: [f64; N],
    c: f64,
    cd: f64,
    cm: f64,
    ip: usize,
    jp: usize,
    gauss: Option<f64>,
}

impl Marsaglia {
    //C*** FIRST CALL RSTART (I, J, K, L)
    //C*** WITH I, J, K, L INTEGERS
    //C*** FROM 1 TO 168, NOT ALL 1
    //C*** NOTE: RSTART CHANGES I, J, K, L
    //C*** SO BE CAREFUL IF YOU REUSE
    //C* * * THEM IN THE CALLING PROGRAM.

    pub fn new(mut i: i32, mut j: i32, mut k: i32, mut l: i32) -> Self {
        let mut u = [0.0; N];
        for e in u.iter_mut().take(N) {
            let mut s = 0.0;
            let mut t = 0.5;
            for _ in 0..24 {
                let m = (((i * j) % 179) * k) % 179;
                i = j;
                j = k;
                k = m;
                l = (53 * l + 1) % 169;
                if l * m % 64 >= 32 {
                    s += t;
                }
                t *= 0.5;
            }
            *e = s;
        }

        Marsaglia {
            u,
            c: 362436. / 16777216.,
            //cd: 7654321. / 167777216.,
            cd: 7654321. / 16777216.,
            cm: 16777213. / 16777216.,
            ip: N - 1,
            jp: 32,
            gauss: None,
        }
    }

    pub fn uni(&mut self) -> f64 {
        let mut uni = self.u[self.ip] - self.u[self.jp];
        if uni < 0.0 {
            uni += 1.0;
        }
        self.u[self.ip] = uni;
        self.ip = if self.ip == 0 { N - 1 } else { self.ip - 1 };
        self.jp = if self.jp == 0 { N - 1 } else { self.jp - 1 };
        self.c -= self.cd;
        if self.c < 0.0 {
            self.c += self.cm;
        }
        uni -= self.c;
        if uni < 0.0 {
            uni + 1.0
        } else {
            uni
        }
    }

    // normally distributed random number
    pub fn gauss(&mut self) -> f64 {
        if let Some(r) = self.gauss {
            self.gauss = None;
            r
        } else {
            let u1 = self.uni();
            let u2 = self.uni();
            let r = (-2.0 * u1.ln()).sqrt();
            let phi = 2.0 * std::f64::consts::PI * u2;
            self.gauss = Some(r * phi.sin());
            r * phi.cos()
        }
    }

    // cauchy distributed random number
    pub fn cauchy(&mut self) -> f64 {
        (2.0 * PI * self.uni()).tan()
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

#[cfg(test)]
mod tests {
    use crate::Marsaglia;
    #[test]
    fn it_works() {
        let mut rng = Marsaglia::new(12, 34, 56, 78);
        let mut l = Vec::new();
        for ii in 1..=20005 {
            let x = rng.uni();
            if ii > 20000 {
                for i in 1..=7 {
                    let v = (x * 16.0_f64.powi(i)) as i32;
                    let v = v % 16;
                    l.push(v);
                }
            }
        }
        #[rustfmt::skip]
        const R: [i32; 35] = [
            6, 3, 11,  3,  0,  4, 0, 
           13, 8, 15, 11, 11, 14, 0, 
            6, 15, 0,  2,  3, 11, 0, 
            5, 14, 2, 14,  4,  8, 0, 
            7, 15, 7, 10, 12,  2, 0,
        ];
        assert!(l.iter().eq(R.iter()));
    }
}
