use marsaglia_rs::{get_float, get_number, yes};
use marsaglia_rs::f::ftest;

// VARIANCE RATIO TEST (F-TEST): COMPARISON OF TWO ERROR BARS.
fn main() {

    println!("VARIANCE RATIO TEST (F-TEST):");
    println!("INPUT: TWO ERROR BARS (NOT VARIANCES!) AND NUMBER OF");
    println!("       INDEPENDENT DATA THEY ARE BASED ON");
    println!("OUTPUT: LIKELIHOOD THAT THEIR RATIO IS DUE TO CHANCE\n");
    loop {
        let eb1 = get_number::<f64>("ENTER ERROR BAR OF 1. DATA POINT:");
        let ndat1 = get_number::<u32>("ENTER NUMBER OF INDEPENDENT MEASUREMENTS FOR 1. DATA POINT:");
        let eb2 = get_number::<f64>("ENTER ERROR BAR OF 2. DATA POINT:");
        let ndat2 = get_number::<u32>("ENTER NUMBER OF INDEPENDENT MEASUREMENTS FOR 2. DATA POINT:");


        let q = ftest(eb1,eb2,ndat1,ndat2);

        println!("Q = {q}");
        println!("OR EQUIVALENTLY: {:.2}%", q * 100.0);
        if !yes("DO YOU WANT ANOTHER RUN? (Y/N).", "OK.", "GOODBYE.") {
            break;
        }
    }
}
