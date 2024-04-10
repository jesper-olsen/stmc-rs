// EMPIRICAL DISTRIBUTION FUNCTION AND Q-DISTRIBUTION.

use gnuplot::{AxesCommon, Caption, Color, Figure };
use marsaglia_rs::marsaglia::Marsaglia;

fn main() {
    const N: usize = 100;

    for title in ["Uniform", "Gauss", "Uniform2"] {
        let mut rng = Marsaglia::new(12, 34, 56, 78);
        let mut data = [0.0; N];
        for e in data.iter_mut() {
            *e = match title {
                "Uniform" => rng.uni(),
                "Gauss"  => rng.gauss(),
                "Uniform2" => 2.0 * (rng.uni() - 0.5),
                _ => unimplemented!()
            };
        }
        data.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mut vx = Vec::new();
        let mut vf = Vec::new();
        let mut vfq = Vec::new();
        vx.push(data[0]);
        vf.push(0.0);
        vfq.push(0.0);
        for j in 0..data.len() - 1 {
            let f = (j + 1) as f64 / data.len() as f64;
            let fq = f.min(1.0 - f);
            vx.push(data[j]);
            vf.push(f);
            vfq.push(fq);
            vx.push(data[j + 1]);
            vf.push(f);
            vfq.push(fq);
        }
        vx.push(data[data.len() - 1]);
        vf.push(1.0);
        vfq.push(0.0);

        let mut fg = Figure::new();
        fg.set_title( format!("Empirical F and F_q - {}", title).as_str());
        fg.axes2d()
            .set_x_label("x", &[])
            .set_y_label("F", &[])
            .lines(&vx, &vf, &[Caption("F"), Color("red")])
            .lines(&vx, &vfq, &[Caption("Fq"), Color("blue")]);
        fg.show().unwrap();
    }
}
