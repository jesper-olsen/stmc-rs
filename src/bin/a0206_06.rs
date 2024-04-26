use stmc_rs::kolm::kolm2_as;
use stmc_rs::marsaglia::Marsaglia;
use stmc_rs::plot::plot2;

fn main() {
    two_xr_a();
    two_xr_b();
}

fn f1x(y: f64) -> f64 {
    if y <= 0.0 {
        0.0
    } else if y <= 1.0 {
        y
    } else {
        1.0
    }
}

fn f2x(y: f64) -> f64 {
    if y <= 0.0 {
        0.0
    } else if y <= 1.0 {
        0.5 * y.powi(2)
    } else if y <= 2.0 {
        0.5 * y * (4.0 - y) - 1.0
    } else {
        1.0
    }
}

fn two_xr_a() {
    const NDAT: usize = 1000;
    const NRPT: usize = 1000;
    let mut a = [0.0f64; NDAT];
    let mut q = [0.0f64; NRPT];
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    for irpt in 0..NRPT {
        for e in a.iter_mut() {
            *e = rng.uni() + rng.uni();
        }
        a.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for e in a.iter_mut() {
            *e = f2x(*e);
        }
        (_, q[irpt]) = kolm2_as(&a);
        if irpt == 0 {
            println!("irpt,Q: {irpt:6}{:9.4}", q[irpt]);
        }
    }
    q.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut v = Vec::new();
    for irpt in 0..NRPT {
        let x = (irpt + 1) as f64 / (NRPT as f64 + 1.0);
        let mut peaked = q[irpt];
        if peaked > 0.5 {
            peaked = 1.0 - peaked;
        }
        //     println!("{x:14.6}{:14.6}{peaked:14.6}", q[irpt]);
        v.push((x, peaked));
    }
    for irpt in 0..NRPT {
        a[irpt] = f1x(q[irpt]);
    }
    let (_, qq) = kolm2_as(&a);
    println!(" Q for the Q distribution = {qq:9.4}");

    let mut graphs = Vec::new();
    graphs.push((String::new(), v));

    plot2(
        "fig_a0206_06a.png",
        "Part 1", // title
        "x",
        "y",
        graphs,
        0.0..1.0,
        0.0..0.5,
    );
}
fn get_rng(irpt: i32) -> Marsaglia {
    let ij = 1 + 1801;
    let kl = 9373 + 0 + irpt;
    let i = (ij / 177) % 177 + 2; // ! I=0 for IJ=0 and IJ=31329.
    let j = ij % 177 + 2;
    let k = (kl / 169) % 178 + 1; // K=0 for KL=0 and KL=30082.
    let m = kl % 169;
    Marsaglia::new(i, j, k, m)
}

fn two_xr_b() {
    const NDAT: usize = 1000;
    const NRPT: usize = 1000;
    let mut a = [0.0f64; NDAT];
    let mut b = [0.0f64; NDAT];
    let mut q = [0.0f64; NRPT];

    for (irpt, e) in q.iter_mut().enumerate() {
        let mut rng = get_rng(irpt as i32);
        for e in a.iter_mut() {
            *e = rng.uni();
        }

        let mut rng = get_rng(irpt as i32 + 1);
        for (i, e) in b.iter_mut().enumerate() {
            *e = rng.uni();
            a[i] += *e;
        }
        a.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for e in a.iter_mut() {
            *e = f2x(*e);
        }
        (_, *e) = kolm2_as(&a);
        if irpt == 0 {
            println!("irpt,Q: {irpt:6}{:9.4}", *e);
        }
    }
    q.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut v = Vec::new();
    for irpt in 0..NRPT {
        let x = (irpt + 1) as f64 / (NRPT as f64 + 1.0);
        let mut peaked = q[irpt];
        if peaked > 0.5 {
            peaked = 1.0 - peaked;
        }
        //      println!("{x:14.6}{:14.6}{peaked:14.6}", q[irpt]);
        v.push((x, peaked));
    }
    for irpt in 0..NRPT {
        a[irpt] = f1x(q[irpt]);
    }
    let (_, qq) = kolm2_as(&a);
    println!(" Q for the Q distribution = {qq:9.4}");

    let mut graphs = Vec::new();
    graphs.push((String::new(), v));

    plot2(
        "fig_a0206_06b.png",
        "Part 2", // title
        "x",
        "y",
        graphs,
        0.0..1.0,
        0.0..0.5,
    );
}
