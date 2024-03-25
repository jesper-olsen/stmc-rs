use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::steb::{steb0, steb2};

// DISASTER WITH SMALL BIN SIZES. AVERAGES WEIGHTED WITH ONE OVER
// THEIR ESTIMATED VARIANCE. RESULTS FOR MEAN VALUES AND ERROR BARS.

const NDAT: usize = 19200;
const KMAX: usize = 6;
const JDAT: [usize; KMAX] = [2, 3, 4, 8, 16, 32];

fn main() {
    let mut dat0 = [0.0; NDAT];
    let mut dat2 = [0.0; NDAT];
    let mut eb = [0.0; NDAT];
    let mut wght = [0.0; NDAT];
    let mut xm = [0.0; KMAX + 1];
    let mut xe = [0.0; KMAX + 1];
    let mut kdat = [0; KMAX];

    let mut rng = Marsaglia::new(12, 34, 56, 78);

    for e in dat0.iter_mut() {
        *e = rng.uni();
    }

    (xm[0], xe[0]) = steb0(&dat0);

    for k in 0..KMAX {
        for i in (1..=NDAT).step_by(JDAT[k]) {
            let j = 1 + i / JDAT[k];
            let z = i - 1;
            (dat2[j - 1], eb[j - 1]) = steb0(&dat0[z..z + JDAT[k]]);
        }
        wght[0] = -1.0;
        kdat[k] = NDAT / JDAT[k];
        (xm[k + 1], xe[k + 1]) = steb2(&dat2[0..kdat[k]], &eb[0..kdat[k]], &mut wght[0..kdat[k]]);
    }

    print!("\nBinsize:       ");
    for e in JDAT {
        print!(" {e:6} ");
    }

    print!("\nData:    {NDAT} ");
    for e in kdat {
        print!(" {:6} ", e);
    }

    print!("\nMeans: ");
    for e in xm {
        print!(" {e:1.4} ");
    }

    print!("\nErrors:");
    for e in xe {
        print!(" {e:1.4} ");
    }
    println!();
}
