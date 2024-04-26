use stmc_rs::chi2::chi2_df;
use stmc_rs::marsaglia::Marsaglia;
use stmc_rs::plot::plot2;

// CHI2 Test for events from a dice versus the expectation.

const NDAT: usize = 1000;
const NRPT: usize = 10_000;
const NHIST: usize = 6;
fn main() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    let mut x = [0.0; NHIST];
    let mut q = [0.0; NRPT];
    let mut qfalse = [0.0; NRPT];

    for e in x.iter_mut() {
        *e = NDAT as f64 / NHIST as f64;
    }
    for i in 0..NRPT {
        let mut hist = [0.0; NHIST];
        for _ in 0..NDAT {
            let idice = (NHIST as f64 * rng.uni()) as usize;
            hist[idice] += 1.0;
        }
        let mut chi2 = 0.0;
        for j in 0..NHIST {
            chi2 += (hist[j] - x[j]).powi(2) / x[j];
        }
        q[i] = 1.0 - chi2_df(chi2, NHIST - 1);
        qfalse[i] = 1.0 - chi2_df(chi2, NHIST);
        //println!("IRPT,CHI2,Q:{i},{chi2},{}", q[i]);
    }
    q.sort_by(|a, b| a.partial_cmp(b).unwrap());
    qfalse.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut graphs = Vec::new();
    let mut v1 = Vec::new();
    let mut v2 = Vec::new();
    for i in 0..NRPT {
        let mut f = (i + 1) as f64 / (NRPT as f64 + 1.0);
        if f > 0.5 {
            f = 1.0 - f;
        }
        println!("{} {} {}", q[i], qfalse[i], f);
        v1.push((q[i], f));
        v2.push((qfalse[i], f));
    }
    graphs.push((String::new(), v1));
    graphs.push((String::new(), v2));

    //reproduces fig 2.3
    plot2(
        "fig_a0206_01.png",
        "\u{1D6D8}\u{00B2} test for dice",
        "Q",
        "F^p",
        graphs,
        0.0..1.0,
        0.0..0.5,
    );
}
