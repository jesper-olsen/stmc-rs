// Rust implementation of the pseudo random number generator presented in
// TOWARD A UNIVERSAL RANDOM NUMBER GENERATOR
// George MARSAGLIA and Arif ZAMAN
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
        for ii in 0..N {
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
            u[ii] = s;
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
}

fn main() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    for ii in 1..=20005 {
        let x = rng.uni();
        if ii > 20000 {
            for i in 1..=7 {
                let v = (x * 16.0_f64.powi(i)) as i32;
                let v = v % 16;
                print!("{v:8}");
            }
            println!();
        }
    }
}
