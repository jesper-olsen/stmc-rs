use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::plot::plot2;
use marsaglia_rs::steb::{datjack, steb0, stebj0, stebj1};

fn main() {
    jack_eb();
}

// ILLUSTRATION THAT JACKKNIFE ERROR BARS:
// (1) JACKKNIFE ERROR BARS AND STANDARD ERROR BARS
//     AGREE FOR UNBIASED QUANTITIES.
// (2) JACKKNIFE ERROR BAR AND BIAS FOR THE MEAN
//     EXPECTATION VALUE SQUARED.
const NDAT: usize = 320;
fn jack_eb() -> f64 {
    let mut data = [0.0; NDAT];
    let mut datj = [0.0; NDAT];
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    for e in data.iter_mut() {
        *e = rng.uni();
    }
    datjack(&data, &mut datj);
    let (xm, xe) = steb0(&data);
    println!("STEB0:  {xm:.7} +/- {xe:.7}");
    let (xm, _, xe) = stebj0(&datj);
    println!("STEBJ0: {xm:.7} +/- {xe:.7}");
    let xmean = data.iter().sum::<f64>() / data.len() as f64;
    println!("NO BIAS!");
    println!("Now f(\\hat{{x}})=\\hat{{x}}^2:");
    let xmean2 = xmean.powi(2);
    for e in datj.iter_mut() {
        *e = (*e).powi(2);
    }
    let (xm, _, xe) = stebj0(&datj);
    println!("STEBJ0: {xm:.7} +/- {xe:.7}");
    let (xjm, xmm2, _, xm2e) = stebj1(&datj, xmean2);
    println!("STEBJ1: XMEAN2={xmean2:.7} +/-{xe:.7}");
    println!("STEBJ1: XMM2=  {xmm2:.7} +/-{xm2e:.7}");
    let biasm = xmean2 - xmm2;
    println!("EST. BIAS =    {biasm:.7}");
    biasm
}

#[cfg(test)]
mod tests {
    use crate::jack_eb;
    #[test]
    fn jack_eb_test() {
        let bias = jack_eb();
        assert_eq!(bias, 0.0002616591225786946);
    }
}
