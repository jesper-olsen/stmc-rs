use marsaglia_rs::chi2::chi2pdf_xq;

// ERROR BAR OF THE VARIANCE: ASYMPTOTIC LARGE N BEHAVIOUR.
// Table 2.9 p. 82

const KMAX: usize = 14;

fn main() {
    println!("  CONFIDENCE INTERVALS FOR VARIANCE ESTIMATE FROM NDAT DATA");
    println!("\n  NDAT=2**K");
    println!("   NDAT  K  \\\\ q   .025       .150       .500       .850       .975");

    for k in 1..=KMAX {
        let ndat = 2u32.pow(k as u32);
        let nf = (ndat - 1) as usize;
        print!("\n  {ndat:5} {k:2}  ");
        for chi in [0.975, 0.85, 0.5, 0.15, 0.025] {
            let s = 1.0 / chi2pdf_xq(chi, nf);
            print!("{s:11.3}");
        }
    }
    println!();
}
