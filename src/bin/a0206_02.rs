use stmc_rs::chi2::chi2_df;
use stmc_rs::marsaglia::Marsaglia;
use stmc_rs::plot::plot2;
use std::process;

// CHI2 Test for events from one dice versus another.
fn main() {
    const NHIST: usize = 6;
    const NDA1: usize = 1000;
    const NDA2: usize = 1000;
    const NRPT: usize = 10000;
    let mut q = [0.0; NRPT];
    let ratio = NDA1 as f64 / NDA2 as f64;

    let mut rng = Marsaglia::new(12, 34, 56, 78);
    for e in q.iter_mut() {
        let mut hist1 = [0; NHIST];
        let mut hist2 = [0; NHIST];
        for _ in 0..NDA1 {
            let idice = (NHIST as f64 * rng.uni()) as usize;
            hist1[idice] += 1;
        }
        for _ in 0..NDA2 {
            let idice = (NHIST as f64 * rng.uni()) as usize;
            hist2[idice] += 1;
        }
        let mut chi2 = 0.0;
        for ihist in 0..NHIST {
            let sum = hist1[ihist] as f64 + ratio.powi(2) * hist2[ihist] as f64;
            if sum < 0.5 {
                println!("CHI2_TEST: Singular, no event.");
                process::exit(1);
            }
            chi2 += (hist1[ihist] as f64 - ratio * hist2[ihist] as f64).powi(2) / sum;
        }
        *e = 1.0 - chi2_df(chi2, NHIST - 1);
    }

    q.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut graphs = Vec::new();
    let mut v = Vec::new();

    for irpt in 0..NRPT {
        let mut f = (irpt + 1) as f64 / (NRPT + 1) as f64;
        if f > 0.5 {
            f = 1.0 - f;
        }
        //println!("{:11.5} {f:11.5}", q[irpt]);
        v.push((q[irpt], f));
    }
    graphs.push((String::new(), v));

    plot2(
        "fig_a0206_02.png",
        "\u{1D6D8}\u{00B2} test for dice",
        "x",
        "y",
        graphs,
        0.0..1.0,
        0.0..0.5,
    );
}
