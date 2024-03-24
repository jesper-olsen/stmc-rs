use marsaglia_rs::gau::gaudif;
use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::plot::plot;

//C Generate 2*10000 Gaussian random numbers, let the difference
//C between the mean expectation values be 0,1,2, ....
//C Perform each time the Gaussian difference test and histogram
//C the resulting Q-distribution.

fn main() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);

    let mut qdat = [0.0f64; 10000];
    let eb1 = 1.0;
    let eb2 = 1.0;

    let xmin = 0.0;
    let xmax = 1.0;
    let mut ymax = 0.0;
    let mut l = Vec::new();
    for del_mean in 0..=4 {
        let del_mean = del_mean as f64;
        for e in qdat.iter_mut() {
            let x1 = rng.gauss();
            let x2 = rng.gauss();
            let x2 = x2 + del_mean;
            *e = gaudif(x1, eb1, x2, eb2);
        }

        let mut hist = [0.0; 10];
        to_hist(true, &mut hist, &qdat, xmin, xmax);

        let (v, ym) = hist2graph(&hist, xmin, xmax);
        ymax = f64::max(ymax, ym);

        l.push((format!("{del_mean}"), v));
    }

    plot("fig.png", "Histogram", "x", "f", l, xmin, xmax, ymax);
}

fn to_hist(normalise: bool, hist: &mut [f64], data: &[f64], xmin: f64, xmax: f64) {
    let factor = hist.len() as f64 / (xmax - xmin);
    let mut n = 0;
    for &x in data {
        if x >= xmin && x <= xmax {
            let i = (factor * (x - xmin)) as usize;
            hist[i] += 1.0;
            n += 1;
        }
    }
    if n != data.len() {
        println!(
            "to_hist: WARNING! Omitted {} of {} samples",
            data.len() - n,
            data.len()
        );
    }

    if normalise && n > 0 {
        let delta = (xmax - xmin) / (hist.len() as f64);
        for e in hist.iter_mut() {
            *e *= 1.0 / n as f64 / delta;
        }
    }
}

fn hist2graph(hist: &[f64], xmin: f64, xmax: f64) -> (Vec<(f64, f64)>, f64) {
    let mut ymax = 0.0;
    let dx = (xmax - xmin) / hist.len() as f64;
    let mut v = Vec::new();
    for (bin, y) in hist.iter().enumerate() {
        let x = bin as f64 / 10.0;
        if *y > ymax {
            ymax = *y
        }
        v.push((x, *y));
        v.push((x + dx, *y));
    }
    //let ymax= hist.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (v, ymax)
}
