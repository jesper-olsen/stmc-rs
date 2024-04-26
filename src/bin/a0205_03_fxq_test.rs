use stmc_rs::f::{f_df, f_xq};

// F-Ratio: test of the correctness of the q-tile function F_XQ.

fn main() {
    let n = 28;
    let ndvd = 7;
    for i in 1..=n {
        let f = i as f64 / ndvd as f64;

        print!("{f:10.5} ");
        let mut afdf = [0.0; 3];
        let mut afq = [0.0; 3];
        for (j, (nf1, nf2)) in [(16, 16), (32, 32), (32, 16)].iter().enumerate() {
            afdf[j] = f_df(f, *nf1, *nf2);
            afq[j] = f_xq(afdf[j], *nf1, *nf2);
        }
        for e in afdf {
            print!("{e:10.5} ");
        }
        for e in afq {
            print!("{e:10.5} ");
        }
        println!();
    }
}
