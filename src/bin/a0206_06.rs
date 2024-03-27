use marsaglia_rs::kolm::kolm2_as;
use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::plot::plot2;
use std::process;

// CHI2 Test for events from one dice versus another.
fn main() {
    two_xr_a();
    //    let mut rng = Marsaglia::new(12, 34, 56, 78);
    //    q.sort_by(|a, b| a.partial_cmp(b).unwrap());
    //
    //    let mut graphs = Vec::new();
    //    let mut v = Vec::new();
    //
    //    for irpt in 0..NRPT {
    //        let mut f = (irpt + 1) as f64 / (NRPT + 1) as f64;
    //        if f > 0.5 {
    //            f = 1.0 - f;
    //        }
    //        //println!("{:11.5} {f:11.5}", q[irpt]);
    //        v.push((q[irpt], f));
    //    }
    //    graphs.push((String::new(), v));
    //
    //    plot2(
    //        "fig_a0206_02.png",
    //        "\u{1D6D8}\u{00B2} test for dice",
    //        "x",
    //        "y",
    //        graphs,
    //        0.0,
    //        1.0,
    //        0.5,
    //    );
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
    const ndat: usize = 1000;
    const nrpt: usize = 1000;
    let mut a = [0.0f64; ndat];
    let mut q = [0.0f64; nrpt];
    let mut rng = Marsaglia::new(12, 34, 56, 78);

    for irpt in 0..nrpt {
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

    for irpt in 0..nrpt {
        let x = (irpt + 1) as f64 / (nrpt as f64 + 1.0);
        let mut peaked = q[irpt];
        if peaked > 0.5 {
            peaked = 1.0 - peaked;
        }
        println!("{x:14.6}{:14.6}{peaked:14.6}", q[irpt]);
    }
    for irpt in 0..nrpt {
        a[irpt] = f1x(q[irpt]);
    }
    let (_, qq) = kolm2_as(&a);
    println!(" Q for the Q distribution = {qq:9.4}");
}
