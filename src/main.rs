// Rust implementation of the pseudo random number generator presented in
// TOWARD A UNIVERSAL RANDOM NUMBER GENERATOR
// George MARSAGLIA, Arif ZAMAN and Wai Wan TSANG
// Statistics & Probability Letters 8 (1990) 35-39

const N: usize = 97;

struct Marsaglia {
    u: [f64; N],
    c: f64,
    cd: f64,
    cm: f64,
    ip: usize,
    jp: usize,
}

impl Marsaglia {
    //C*** FIRST CALL RSTART (I, J, K, L)
    //C*** WITH I, J, K, L INTEGERS
    //C*** FROM 1 TO 168, NOT ALL 1
    //C*** NOTE: RSTART CHANGES I, J, K, L
    //C*** SO BE CAREFUL IF YOU REUSE
    //C* * * THEM IN THE CALLING PROGRAM.
    fn new(mut i: i32, mut j: i32, mut k: i32, mut l: i32) -> Self {
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
        }
    }

    fn uni(&mut self) -> f64 {
        let mut uni = self.u[self.ip] - self.u[self.jp];
        if uni < 0.0 {
            uni += 1.0;
        }
        self.u[self.ip] = uni;
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
}

fn main() {
    let p = 0.4;
    println!("\nFrequency of samples in [0;{p})");
    println!("k, #samples, frequency, error");
    for k in 1..12 {
        let m = 2_u64.pow(2 * k - 1);
        let mut rng = Marsaglia::new(12, 34, 56, 78);
        let n = (0..m).filter(|_| rng.uni() < p).count();
        let r = (n as f64) / (m as f64);
        println!("{k:3} {m:8} {r:.20} {:+e}", (r - p).abs());
    }

    use std::collections::HashMap;
    for m in [100, 100000] {
        let mut rng = Marsaglia::new(12, 34, 56, 78);
        let mut histogram: HashMap<usize, usize> = HashMap::new();

        (0..m).map(|_| (10.0 * rng.uni()) as usize).for_each(|bin| {
            let count = histogram.entry(bin).or_insert(0);
            *count += 1
        });

        println!("\nHistogram - {m} samples, 10 bins");
        let mut bins: Vec<usize> = histogram.keys().cloned().collect();
        bins.sort();
        let width = (m as f64).ln() as usize;
        for bin in bins {
            println!("Bin {bin}: {:width$}", histogram[&bin]);
        }
    }

    let mut rng = Marsaglia::new(12, 34, 56, 78);
    let mut f1 = rng.uni();
    let mut c1 = 1;
    let mut f2 = f1;
    let mut c2 = 1;
    const M: u64 = 10_000_000_000;
    println!("\nmin, max and their frequencies over {M} samples");

    for _ in 0..M - 1 {
        let f = rng.uni();
        if f == f1 {
            c1 += 1;
        } else if f < f1 {
            f1 = f;
            c1 = 0;
        }
        if f == f2 {
            c2 += 1;
        } else if f > f2 {
            f2 = f;
            c2 = 0;
        }
    }
    println!("min {f1}; count {c1}");
    println!("max {f2}; count {c2}");
}

#[cfg(test)]
mod tests {
    #[cfg(test)]
    use crate::Marsaglia;
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
