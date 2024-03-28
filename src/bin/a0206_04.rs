use marsaglia_rs::cau::cau_df;
use marsaglia_rs::gau::gau_df;
use marsaglia_rs::kolm::{kolm1, kolm1_as};
use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::plot::plot2;

// Statistical investigation of the 1-sided Kolmogorov test.
// Examples: Uniform, Gaussian and Cauchy random numbers.f

fn main() {
    a020604(true);
    a020604(false);
}

fn a020604(part1: bool) {
    //const N: usize = 5;
    const N: usize = 1000;
    const NRPT: usize = 10_000;
    let mut data = [0.0; N];
    let mut fxct = [0.0; N];
    let mut q1 = [0.0f64; NRPT];
    let mut q2 = [0.0f64; NRPT];
    let mut q3 = [0.0f64; NRPT];
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    let distribution = 1;

    println!(
        "{} random numbers.",
        match distribution {
            1 => "Cauchy",
            2 => "Gauss",
            _ => "Uniform",
        }
    );

    for irpt in 0..NRPT {
        for e in data.iter_mut() {
            *e = match distribution {
                1 => rng.cauchy(),
                2 => rng.gauss(),
                _ => rng.uni(),
            }
        }

        if N > 1 {
            data.sort_by(|a, b| a.partial_cmp(b).unwrap());
        }
        for i in 0..data.len() {
            fxct[i] = match distribution {
                1 => cau_df(data[i]),
                2 => gau_df(data[i]),
                _ => data[i], // uniform
            }
        }
        if part1 {
            (_, _, q1[irpt], q2[irpt]) = kolm1(&fxct);
        } else {
            (_, _, q1[irpt], q2[irpt]) = kolm1_as(&fxct);
        }
        q3[irpt] = q1[irpt].min(q2[irpt]);
    }

    q1.sort_by(|a, b| a.partial_cmp(b).unwrap());
    q2.sort_by(|a, b| a.partial_cmp(b).unwrap());
    q3.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut graphs = Vec::new();
    let mut v1 = Vec::new();
    let mut v2 = Vec::new();
    let mut v3 = Vec::new();

    for irpt in 0..NRPT {
        let mut f = (irpt + 1) as f64 / (NRPT + 1) as f64;
        if f > 0.5 {
            f = 1.0 - f;
        }
        println!(
            "{:11.5}{:11.5}{:11.5}{f:11.5}",
            q1[irpt], q2[irpt], q3[irpt]
        );
        v1.push((q1[irpt], f));
        v2.push((q2[irpt], f));
        v3.push((q3[irpt], f));
    }

    graphs.push((String::new(), v1));
    graphs.push((String::new(), v2));
    graphs.push((String::new(), v3));

    plot2(
        if part1 {
            "fig_a0206_04a.png"
        } else {
            "fig_a0206_04b.png"
        },
        if part1 { "kolm1" } else { "kolm1_as" }, // title
        "x",
        "y",
        graphs,
        0.0,
        1.0,
        0.5,
    );
}
