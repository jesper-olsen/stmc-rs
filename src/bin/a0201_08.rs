use marsaglia_rs::gaudif;
use marsaglia_rs::marsaglia::Marsaglia;
use plotters::prelude::*;

//C Generate 2*10000 Gaussian random numbers, let the difference
//C between the mean expectation values be 0,1,2, ....
//C Perform each time the Gaussian difference test and histogram
//C the resulting Q-distribution.

fn main() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);

    const NHIST: usize = 10;
    //let del_mean = 0.0;
    let del_mean = 4.0;
    let mut qdat = [0.0f64; 10000];
    let eb1 = 1.0;
    let eb2 = 1.0;

    let xmin = 0.0;
    let xmax = 1.0;
    let mut ymax = 0.0;
    let mut l = Vec::new();
    for del_mean in 0..4 {
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
        ymax = f64::max(ymax,ym);

        l.push((format!("{del_mean}"), v));
    }

    plot("fig.png", "Histogram", l, xmin, xmax, ymax);
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

    if normalise {
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
    for bin in 0..hist.len() {
        let x = bin as f64 / 10.0;
        let y = hist[bin];
        if y > ymax {
            ymax = hist[bin];
        }
        v.push((x, y));
        v.push((x + dx, y));
    }
    //let ymax= hist.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (v, ymax)
}

fn plot(
    fname: &str,
    caption: &str,
    graphs: Vec<(String, Vec<(f64, f64)>)>,
    xmin: f64,
    xmax: f64,
    ymax: f64,
) {
    println!("Saving {fname}");
    let root_area = BitMapBackend::new(fname, (600, 400)).into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_area)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .caption(caption, ("sans-serif", 40))
        //.build_cartesian_2d(-3f64..3f64, 0f64..ymax)
        .build_cartesian_2d(xmin..xmax, 0f64..ymax)
        .unwrap();

    ctx.configure_mesh().x_desc("x").y_desc("f").draw().unwrap();

    graphs.into_iter().enumerate().for_each(|(i, (label, h))| {
        let colour = match i {
            0 => RED,
            _ => BLUE,
        };
        ctx.draw_series(
            AreaSeries::new(
                h,
                0.0,             // Baseline
                colour.mix(0.2), // Make the series opac
            )
            .border_style(colour), // Make a brighter border
        )
        .unwrap()
        .label(label)
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], colour));
    });

    ctx.configure_series_labels()
        .border_style(BLACK)
        .background_style(WHITE.mix(0.8))
        .position(SeriesLabelPosition::UpperRight)
        .margin(20)
        .draw()
        .unwrap();
}
