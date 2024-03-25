use marsaglia_rs::gau::{gau_pd, gau_qdf};
use marsaglia_rs::plot::plot2;
use marsaglia_rs::student::{stud_pd, stud_qdf};

const NGNU: i32 = 60;
const NFMAX: usize = 8;
const NFILES: usize = 2;
const MF: usize = NFMAX / NFILES;

/// PROBABILITY DENSITIES AND Q-DISTRIBUTION FUNCTION
/// FOR CHI2 PER DEGREE OF FREEDOM (PDF) DISTRIBUTION.
fn main() {
    let mut graphs1 = Vec::new(); // densities
    let mut graphs2 = Vec::new(); // Q distribution
    for i in 0..NFILES {
        for j in 0..MF {
            let nf = 1 + j + i * MF;
            let mut v1 = Vec::new();
            let mut v2 = Vec::new();
            let mut v3 = Vec::new();
            let mut v4 = Vec::new();
            for z in -NGNU..=NGNU {
                let t = 3.0 * z as f64 / NGNU as f64;
                v1.push((t, stud_pd(t, nf)));
                v3.push((t, stud_qdf(t, nf)));
                if i == NFILES - 1 && j == MF - 1 {
                    v2.push((t, gau_pd(t, 0.0, 1.0)));
                    v4.push((t, gau_qdf(t)));
                }
            }

            graphs1.push((String::new(), v1));
            if !v2.is_empty() {
                graphs1.push((String::new(), v2))
            }
            graphs2.push((String::new(), v3));
            if !v4.is_empty() {
                graphs2.push((String::new(), v4))
            }
        }
    }

    plot2(
        "fig_a0203_02a.png",
        "Normal and Student PDFs (1-8 dg of freedom)",
        "t",
        "f",
        graphs1,
        -3.0,
        3.0,
        0.4,
    );
    plot2(
        "fig_a0203_02b.png",
        "Normal and student Q dfs (1-8 dg of freedom)",
        "t",
        "fq",
        graphs2,
        -3.0,
        3.0,
        0.5,
    );
}
