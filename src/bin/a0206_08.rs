use marsaglia_rs::kolm::kolm2_as2;
use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::plot::plot2;

fn main() {
    k2();
}

// Kolmogorov DEL-distribution for two empirical uniform distributions.
fn k2() {
    let mut graphs = vec![];
    for (n1, n2, nrpt) in [
        (72, 96, 1_000_000),
        //(12, 16, 1_000_000),
        //(12, 16, 100),
    ] {
        let mut dat1 = vec![0.0f64; n1];
        let mut dat2 = vec![0.0f64; n2];
        let mut del = vec![0.0f64; nrpt + 1];
        let mut qas = vec![0.0f64; nrpt + 2];
        let distribution = 0;
        println!(
            "{} random numbers.",
            match distribution {
                1 => "Cauchy",
                2 => "Gauss",
                _ => "Uniform",
            }
        );
        println!("{n1} {n2}");
        let mut rng = Marsaglia::new(12, 34, 56, 78);
        for irpt in 1..nrpt + 1 {
            for e in dat1.iter_mut().chain(dat2.iter_mut()) {
                *e = match distribution {
                    1 => rng.cauchy(),
                    2 => rng.gauss(),
                    _ => rng.uni(),
                }
            }
            dat1.sort_by(|a, b| a.partial_cmp(b).unwrap());
            dat2.sort_by(|a, b| a.partial_cmp(b).unwrap());
            (del[irpt], qas[irpt]) = kolm2_as2(&dat1, &dat2);
        }
        del[1..].sort_by(|a, b| a.partial_cmp(b).unwrap());
        qas[1..=nrpt].sort_by(|a, b| a.partial_cmp(b).unwrap());

        qas[nrpt + 1] = 1.0;
        let mut q0 = 1.0;
        let mut q = 0.0;
        println!("{0:10}{:12.6}{q0:12.6}{:12.6}", del[0], qas[nrpt + 1]);
        let mut v = vec![(del[0], q0)];
        for irpt in 2..=nrpt {
            if del[irpt] > del[irpt - 1] {
                let f = irpt as f64 / (nrpt + 1) as f64;
                q = 1.0 - f;
                let accuracy = (qas[nrpt + 2 - irpt] - q0).abs() / q0;
                println!(
                    "{irpt:10}{:12.6}{q0:12.6}{:12.6}{accuracy:12.6}",
                    del[irpt - 1],
                    qas[nrpt + 2 - irpt]
                );
                let accuracy = (qas[nrpt + 1 - irpt] - q).abs() / q;
                println!(
                    "{irpt:10}{:12.6}{q:12.6}{:12.6}{accuracy:12.6}",
                    del[irpt - 1],
                    qas[nrpt + 1 - irpt]
                );
                v.extend([(del[irpt - 1], q0), (del[irpt - 1], q)]);
                q0 = q;
            }
        }
        println!("{nrpt:10}{:12.6}{q0:12.6}{:12.6}", del[nrpt], qas[1]);
        println!("{nrpt:10}{:12.6}{q:12.6}{:12.6}", del[nrpt], qas[1]);
        v.extend([(del[nrpt], q0), (del[nrpt], q)]);
        graphs.push((String::new(), v));
    }
    plot2(
        "fig_a0206_08.png",
        "", // title
        "Q",
        "delta",
        graphs,
        0.0..1.0,
        0.0..1.0,
    );
}
