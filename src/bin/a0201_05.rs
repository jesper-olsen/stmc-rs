use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::steb::{steb0, steb2};

fn main() {
    const KMAX: usize = 10;
    const NMAX: usize = 2usize.pow(KMAX as u32 - 1);
    let mut data = [0.0f64; KMAX];
    let mut dat1 = [0.0f64; NMAX];
    let mut eb = [0.0f64; KMAX];
    let mut wght = [0.0f64; KMAX];
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    wght[0] = -1.0;
    for k in 0..KMAX {
        let ndat = if k == KMAX - 1 {
            2
        } else {
            2usize.pow(1 + k as u32)
        };
        let a = &mut dat1[0..ndat];
        for e in a.iter_mut() {
            *e = rng.uni();
        }
        (data[k], eb[k]) = steb0(a);
    }
    let (xm, xe) = steb2(&data, &eb, &mut wght);
    println!("UNIFORM RANDOM NUMBERS: XM={xm:.7} +/- {xe:.7}");
}
