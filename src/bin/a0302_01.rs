use gnuplot::{Caption, Color, Figure};
use marsaglia_rs::gau::{gau_df, gau_qdf};
use marsaglia_rs::marsaglia::Marsaglia;

use std::fs::File;
use std::io::prelude::*;

// Metropolis Generation of Gaussian  Random numbers.

const NDAT: usize = 2000;
//const NAT: usize = 20000;

fn main() {
    const A: f64 = 3.0;
    let mut data = [0.0f64; NDAT];
    println!("Metropolis Generation of Gaussian random numbers");
    println!(" NDAT,A:               {NDAT:10} {A:12.2}");
    let acpt = gau_metro(A, &mut data);
    println!(" Acceptance rate found {acpt:22.6}");

    data.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // df_gnu
    let mut adf = [vec![data[0]], vec![0.0], vec![0.0]];
    for i in 0..data.len() - 1 {
        let f = (i + 1) as f64 / data.len() as f64;
        let fq = f.min(1.0 - f);
        adf[0].extend(&data[i..i + 2]);
        adf[1].extend(&[f, f]);
        adf[2].extend(&[fq, fq]);
    }
    adf[0].push(*data.last().unwrap());
    adf[1].push(1.0);
    adf[2].push(0.0);

    let mut file = File::create("df01.d").unwrap();
    for ((x, f), fq) in adf[0].iter().zip(adf[1].iter()).zip(adf[2].iter()) {
        writeln!(file, "{x:16.7}{f:16.7}{fq:16.7}").unwrap();
    }

    let mut fg = Figure::new();
    fg.set_title("df01");
    fg.axes2d()
        .lines(&adf[0], &adf[1], &[Caption("Line 1"), Color("red")])
        .lines(&adf[0], &adf[2], &[Caption("Line 2"), Color("blue")]);
    fg.show().unwrap();

    let mut adf2 = [vec![data[0]], vec![0.0], vec![0.0]];
    const NGNU: i32 = 500;
    for ignu in -NGNU..=NGNU {
        let x = 3.0 * ignu as f64 / NGNU as f64;
        adf2[0].push(x);
        adf2[1].push(gau_df(x));
        adf2[2].push(gau_qdf(x));
    }

    let mut file = File::create("dfexact.d").unwrap();
    for ((x, f), fq) in adf[0].iter().zip(adf[1].iter()).zip(adf[2].iter()) {
        writeln!(file, "{x:16.7}{f:16.7}{fq:16.7}").unwrap();
    }

    let mut fg = Figure::new();
    fg.set_title("df peaked");
    fg.axes2d()
        .lines(&adf[0], &adf[2], &[Caption("Line 1"), Color("red")])
        .lines(&adf2[0], &adf2[2], &[Caption("Line 2"), Color("blue")]);
    fg.show().unwrap();

    let mut fg = Figure::new();
    fg.set_title("df all");
    fg.axes2d()
        .lines(&adf[0], &adf[1], &[Caption("Line 1"), Color("red")])
        .lines(&adf[0], &adf[2], &[Caption("Line 2"), Color("green")])
        .lines(&adf2[0], &adf2[1], &[Caption("Line 3"), Color("blue")])
        .lines(&adf2[0], &adf2[2], &[Caption("Line 4"), Color("pink")]);
    fg.show().unwrap();
}

/// METROPOLIS GENERATION OF GAUSSIAN RANDOM NUMBERS.
/// a:            STEPSIZE PARAMETER; DEL \IN [-a,+a].
fn gau_metro(a: f64, data: &mut [f64]) -> f64 {
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    let mut acpt = 0.0;
    let mut x = 0.0f64;
    for e in data.iter_mut() {
        let x2 = x.powi(2);
        let xp = x + 2.0 * a * (rng.uni() - 0.5);
        let xp2 = xp.powi(2);
        if xp2 <= x2 || (-0.5 * (xp2 - x2)).exp() >= rng.uni() {
            x = xp;
            acpt += 1.0;
        }
        *e = x;
    }
    acpt / data.len() as f64
}
