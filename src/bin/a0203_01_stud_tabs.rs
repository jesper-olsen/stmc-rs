use marsaglia_rs::gauss_cdf;
use marsaglia_rs::student::{stud_df, stud_xq};

///  STUDENT DISTRIBUTION: CONFIDENCE LEVELS FOR SMALL GAUSSIAN SAMPLES.

const NSI: usize = 5;
const MDAT: usize = 64;

fn main() {
    println!("\n   Table 1)\n");
    print!("    N \\ S ");
    for i in 1..=NSI {
        print!("{:11.5} ", i as f64);
    }
    println!("\n");

    for i in 2..=MDAT {
        let nf = i - 1;
        print!("{i:9} ");
        for j in 1..=NSI {
            let t = -1.0 * j as f64;
            let p = 1.0 - 2.0 * stud_df(t, nf);
            print!("{p:11.5} ");
        }
        println!();
    }

    print!("\n INFINITY ");
    for j in 1..=NSI {
        let t = -1.0 * j as f64;
        let p = 1.0 - 2.0 * gauss_cdf(t);
        print!("{p:11.5} ");
    }
    println!();

    println!("\n   Table 2)\n");
    print!("    N \\ P ");
    for i in 1..=NSI {
        let p = 1.0 - 2.0 * gauss_cdf(-1.0 * i as f64);
        print!("{:11.5} ", p);
    }

    println!("\n");

    for i in 2..=MDAT {
        let nf = i - 1;
        let nj = if i == 2 { 4 } else { NSI };

        print!("{i:9} ");
        for j in 1..=nj {
            let xj = j as f64;
            let p = 1.0 - gauss_cdf(-xj);
            let s = stud_xq(p, nf);
            print!("{s:11.5} ");
        }
        println!();
    }

    print!("\n INFINITY ");
    for j in 1..=NSI {
        let p = j as f64;
        print!("{p:11.5} ");
    }
    println!();
}
