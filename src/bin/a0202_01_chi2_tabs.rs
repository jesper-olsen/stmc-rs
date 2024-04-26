use stmc_rs::chi2::{chi2_df, chi2_xq, chi2pdf_df, chi2pdf_xq};

const NMAX: usize = 8;

fn main() {
    let q = [0.025, 0.15, 0.50, 0.85, 0.975];

    println!("\n CHI2 distribution function\n");
    print!("   CHI2 \\ N ");
    for n in 1..=NMAX {
        print!(" {n:4}  ");
    }
    println!();
    for i in 1..=NMAX {
        let chi2 = i as f64;
        print!("   {chi2:1.2}    ");
        for n in 0..NMAX {
            let f = chi2_df(chi2, n + 1);
            print!(" {:2.3} ", f);
        }
        println!();
    }

    println!("\n Fractiles for the CHI2 distribution\n");
    print!("   N \\ q ");
    for qn in q.iter() {
        print!(" {: >8.3}  ", qn);
    }
    println!("\n");
    for n in 1..=NMAX {
        print!("   {n}     ");
        for e in &q {
            print!(" {: >8.3}  ", chi2_xq(*e, n));
        }
        println!();
    }

    println!("\n CHI2 per degree of freedom distribution function");
    print!("   CHI2 \\ N ");
    for n in 1..=NMAX {
        print!(" {n:4}  ");
    }
    println!();
    for i in 1..=NMAX {
        let chi2 = i as f64 / 4.0;
        print!("   {chi2:1.2}    ");
        for n in 0..NMAX {
            print!(" {:2.3} ", chi2pdf_df(chi2, n + 1));
        }
        println!();
    }

    println!("\n Fractiles for the CHI2 per degree of freedom distribution\n");
    print!("   N \\ q ");
    for qn in &q {
        print!(" {: >8.3}  ", qn);
    }
    println!("\n");
    for n in 1..=NMAX {
        print!("   {n}     ");
        for v in &q {
            print!(" {: >8.3}  ", chi2pdf_xq(*v, n));
        }
        println!();
    }
}
