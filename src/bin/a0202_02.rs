use marsaglia_rs::chi2::{chi2pdf_pd, chi2pdf_qdf};
use marsaglia_rs::plot::plot2;

const NGNU: usize = 60;
const NFMAX: usize = 20;
const NFILES: usize = 5;
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
            for z in 0..NGNU {
                let chi2 = 2.0 * (z + 1) as f64 / NGNU as f64;
                let y = chi2pdf_pd(chi2, nf);
                v1.push((chi2, y));
                let y = chi2pdf_qdf(chi2, nf);
                v2.push((chi2, y));
                //println!("{chi2},{y}");
            }
            graphs1.push((String::new(), v1));
            graphs2.push((String::new(), v2))
        }
    }
    plot2(
        "fig2_1.png",
        "PDFs with 1-20 degrees of freedom",
        "chi square",
        "f",
        graphs1,
        0.0..2.0,
        0.0..1.4,
    );
    plot2(
        "fig2_2.png",
        "Peaked distribution functions with 1-20 degrees of freedom",
        "chi square",
        "fq",
        graphs2,
        0.0..2.0,
        0.0..0.55,
    );
}
