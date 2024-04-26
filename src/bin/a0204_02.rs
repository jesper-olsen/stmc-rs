use stmc_rs::gau::{sebar_e, sebar_e_as};

///ERROR BAR OF THE VARIANCE: ASYMPTOTIC LARGE N BEHAVIOUR.

fn main() {
    const NMAX: usize = 14;
    const PC: f64 = 0.95;

    let mut ndat = 1;
    println!("    NDAT   LOWER VARIANCE      ASYMPTOTIC  UPPER VARIANCE      ASYMPTOTIC");
    for _ in 1..=NMAX {
        ndat *= 2;
        let (ebup1, ebdo1) = sebar_e(ndat, PC);
        let (ebup2, ebdo2) = sebar_e_as(ndat, PC);
        println!(
            "{ndat:8}, {:15.4} {:15.4} {:15.4} {:15.4}",
            ebdo1.powi(2),
            ebdo2.powi(2),
            ebup1.powi(2),
            ebup2.powi(2)
        );
    }
}
