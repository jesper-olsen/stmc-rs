use stmc_rs::gau::gaudif;
use stmc_rs::get_number;

fn main() {
    println!("DIFFERENCE TEST: COMPARISON OF TWO MEANS.");
    println!("INPUT: TWO MEAN VALUES AND THEIR ERROR BARS (NOT VARIANCES!");
    println!("LIKELIHOOD FOR THE OBSERVED DISCREPANCE TO BE DUE TO CHANCE");
    let xm1 = get_number::<f64>("ENTER 1. MEAN VALUE (DATA POINT):");
    let eb1 = get_number::<f64>("ENTER ERROR BAR OF 1. DATA POINT:");
    let xm2 = get_number::<f64>("ENTER 2. MEAN VALUE (DATA POINT):");
    let eb2 = get_number::<f64>("ENTER ERROR BAR OF 2. DATA POINT:");

    let q = gaudif(xm1, eb1, xm2, eb2);

    println!("Q = {q}");
    println!("OR EQUIVALENTLY: {:.2}%", q * 100.0);
}
