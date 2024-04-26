use stmc_rs::marsaglia::Marsaglia;
use stmc_rs::steb::steb0;

fn main() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);

    let mut data = [0.0f64; 1024];
    for e in data.iter_mut() {
        *e = rng.uni();
    }
    let (xm, xe) = steb0(&data);
    println!("UNIFORM  RANDOM NUMBERS: {xm} {xe}");

    let mut data = [0.0f64; 2048];
    for e in data.iter_mut() {
        *e = rng.gauss();
    }
    let (xm, xe) = steb0(&data);
    println!("GAUSSIAN RANDOM NUMBERS: {xm} {xe}");

    let mut data = [0.0f64; 1024];
    for e in data.iter_mut() {
        *e = rng.cauchy();
    }
    let (xm, xe) = steb0(&data);
    println!("CAUCHY   RANDOM NUMBERS: {xm} {xe}");
}
