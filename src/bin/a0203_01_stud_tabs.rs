use marsaglia_rs::gauss_pdf;
use marsaglia_rs::student::{stud_df, stud_xq};

///  STUDENT DISTRIBUTION: CONFIDENCE LEVELS FOR SMALL GAUSSIAN SAMPLES.

const NSI: usize = 5;
const MDAT: usize = 64;

fn main() {
    let mut nf: usize;
    let mut p: [f64; NSI] = [0.0; NSI];
    let mut s: [f64; NSI] = [0.0; NSI];

    // Confidence Level
    for i in 0..NSI {
        s[i] = (i + 1) as f64;
    }

    println!();
    println!("   Table 1)");
    println!();

    println!(
        "{:3}  \\ S {:11.5} {:11.5} {:11.5} {:11.5} {:11.5}",
        "N", s[0], s[1], s[2], s[3], s[4]
    );
    println!();

    for i in 2..=MDAT {
        nf = i - 1;
        for j in 0..NSI {
            p[j] = 1.0 - 2.0 * stud_df(-s[j], nf);
        }
        println!(
            "{:5} {:11.5} {:11.5} {:11.5} {:11.5} {:11.5}",
            i, p[0], p[1], p[2], p[3], p[4]
        );
    }
}
