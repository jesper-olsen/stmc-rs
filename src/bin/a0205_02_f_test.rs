use stmc_rs::f::ftest;
use stmc_rs::{get_number, yes};

// VARIANCE RATIO TEST (F-TEST): COMPARISON OF TWO ERROR BARS.
fn main() {
    println!("VARIANCE RATIO TEST (F-TEST):");
    println!("INPUT: TWO ERROR BARS (NOT VARIANCES!) AND NUMBER OF");
    println!("       INDEPENDENT DATA THEY ARE BASED ON");
    println!("OUTPUT: LIKELIHOOD THAT THEIR RATIO IS DUE TO CHANCE\n");
    loop {
        let eb1 = get_number::<f64>("ENTER ERROR BAR OF 1. DATA POINT:");
        let ndat1 =
            get_number::<u32>("ENTER NUMBER OF INDEPENDENT MEASUREMENTS FOR 1. DATA POINT:");
        let eb2 = get_number::<f64>("ENTER ERROR BAR OF 2. DATA POINT:");
        let ndat2 =
            get_number::<u32>("ENTER NUMBER OF INDEPENDENT MEASUREMENTS FOR 2. DATA POINT:");

        let q = ftest(eb1, eb2, ndat1, ndat2);

        println!("Q = {q}");
        println!("OR EQUIVALENTLY: {:.2}%", q * 100.0);
        if !yes("DO YOU WANT ANOTHER RUN? (Y/N).", "OK.", "GOODBYE.") {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use stmc_rs::f::ftest;
    #[test]
    fn f_test() {
        // Table 2.11, pp.88
        for (eb1, ndat1, eb2, ndat2, q) in [
            (1.0, 16, 1.0, 16, 0.9999999998352209),
            (1.0, 16, 1.0, 8, 0.35927326002873716),
            (1.0, 16, 1.0, 4, 0.27983973897718917),
            (1.0, 32, 1.0, 8, 0.06328268138077053),
            (1.0, 64, 1.0, 16, 0.004584323015007774),
            (1.0, 1024, 1.05, 1024, 0.11879954967334894),
            (1.0, 2048, 1.05, 2048, 0.027331913291198884),
            (1.0, 32, 2.0, 8, 0.8998822630202368),
            (1.0, 1024, 2.0, 256, 0.9841754392334929),
            (1.0, 16, 2.0, 16, 0.010889547060149107),
        ] {
            let qa = ftest(eb1, eb2, ndat1, ndat2);
            println!("{eb1},{ndat1},{eb2},{ndat2}, {q}, {qa}");
            assert_eq!(q, qa);
        }
    }
}
