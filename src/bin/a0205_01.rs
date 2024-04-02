use marsaglia_rs::f::{f_df, f_pd, f_qdf};
use marsaglia_rs::plot::plot2;

const NGNU: usize = 280;
const NDVD: usize = 70;
const NFS: [(u32, u32); 3] = [(16, 16), (32, 32), (32, 16)];

fn main() {
    // F-Ratio distribution function and peaked distribution function.
    let mut graphs = Vec::new(); // densities
    let mut av = [
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ];
    for i in 0..NGNU {
        let f = i as f64 / NDVD as f64;
        for (j, (nf1, nf2)) in NFS.iter().enumerate() {
            av[j * 2].push((f, f_df(f, *nf1, *nf2)));
            av[j * 2 + 1].push((f, f_qdf(f, *nf1, *nf2)));
        }
    }
    for v in av {
        graphs.push((String::new(), v));
    }
    plot2(
        "fig_a0205_01a.png",
        "F-Ratio and peaked distribution functions for (16,16), (32,32) & (32,16) DoF ",
        "F",
        "Fqdf",
        graphs,
        0.0..4.0,
        0.0..1.0,
    );

    // F-Ratio probability densities.
    let mut graphs = Vec::new(); // densities
    let mut av = [Vec::new(), Vec::new(), Vec::new()];
    let mut fsum = 0.0;
    for i in 0..NGNU {
        let f = i as f64 / NDVD as f64;
        for (j, (nf1, nf2)) in NFS.iter().enumerate() {
            let y = f_pd(f, *nf1, *nf2);
            av[j].push((f, y));
            if j == 1 {
                fsum += y / NDVD as f64;
            }
        }
    }

    println!("Fsum={fsum}");
    for v in av {
        graphs.push((String::new(), v));
    }
    plot2(
        "fig_a0205_01b.png",
        "F-Ratio probabilities for (16,16), (32,32) & (32,16) degrees of freedom ",
        "F",
        "Fpd",
        graphs,
        0.0..4.0,
        0.0..1.2,
    );
}
