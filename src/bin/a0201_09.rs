use marsaglia_rs::gau::gaudif;
use marsaglia_rs::steb::steb2;

const Q0CUT: f64 = 0.05;
const P0CUT: f64 = 1.0 - Q0CUT;
const XM0: f64 = 1.75; // Theoretical value.
const EB0: f64 = 0.0;
const LTEST: bool = false;
//const LTEST: bool = true;

fn main() {
    let mut dat = [
        1.98, 1.61, 1.77, 1.79, 1.72, 1.82, 2.24, 2.73, 1.70, 2.01, 1.70, 1.77, 1.85, 1.57, 1.87,
    ];
    let mut eb = [
        0.16, 0.40, 0.40, 0.37, 0.15, 0.20, 0.10, 0.31, 0.43, 0.27, 0.10, 0.20, 0.21, 0.08, 0.30,
    ];
    let mut w = [0.0; 15];

    println!("\nInput data and their error bars:\n");
    for (i, (x, eb)) in dat.iter().zip(eb.iter()).enumerate() {
        println!("{i:6} {x:.2} {eb:.2}", i = i + 1);
    }

    for n in (2..=dat.len()).rev() {
        w[0] = -1.0;
        let (xm, xe) = steb2(&dat[0..n], &eb[0..n], &mut w[0..n]);
        let q = gaudif(XM0, EB0, xm, xe);
        println!("{n:4} {xm:6.3} +/- {xe:6.3} Q={q:6.3}");
        if LTEST {
            println!("\n");
        }
        let nm1 = n - 1;
        let pln = P0CUT.ln() / n as f64;
        let qcut = 1.0 - pln.exp(); // Qcut for nm1 Gaussian difference tests
        let mut qmin = 1.0;
        let mut imin = 0;
        for i in 0..n {
            let xm1 = dat[i];
            let eb1 = eb[i];
            dat[i] = dat[n - 1];
            eb[i] = eb[n - 1];
            w[0] = -1.0;
            let (xm, xe) = steb2(&dat[0..nm1], &eb[0..nm1], &mut w[0..nm1]);
            let q = gaudif(xm1, eb1, xm, xe);
            if LTEST {
                println!("{i:10} {q} {xm1} {eb1:.4} {xm} {xe}");
            }
            dat[i] = xm1;
            eb[i] = eb1;
            if q < qmin {
                qmin = q;
                imin = i;
            }
        }
        println!("\nN,Qcut,Qmin = {n}, {qcut}, {qmin}:");
        if qmin < qcut {
            println!(
                "Imin = {imin} i.e. XM = {} +/- {} eliminated.",
                dat[imin], eb[imin]
            );
            dat[imin] = dat[n - 1];
            eb[imin] = eb[n - 1];
        } else {
            println!(
                "imin={imin} i.e. XM={} +/- {} no data point eliminated.\n",
                dat[imin], eb[imin]
            );
            break;
        }
    }
}
