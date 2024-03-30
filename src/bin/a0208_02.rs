use marsaglia_rs::gamma::gamma_p;
use std::fmt;

fn main() {
    lfit();
}

//LINEAR FIT: FIT OF INPUT DATA BY STRAIGHT LINE Y=A*X+B.
//
//INPUT:   GAUSSIAN DATA          Y (I), I=1,...,N
//         AND THEIR ERROR BARS   YE(I), I=1,...,N.
//
//OUTPUT:  CONSTANTS  A(1), A(2) WITH ERRORS BARS SGA(1), SGA(2).
//         CORRESPONDING TO THE FIT Y=A(1)+A(2)*X. IF N.GE.3:
//         ALSO  CHI2  AND  THE GOODNESS OF FIT Q ARE RETURNED.

// Gaussian data: x, y, error bars
static DATA: [(f64, f64, f64); 6] = [
    (16.0, 0.1086, 0.0007),
    (24.0, 0.1058, 0.0008),
    (34.0, 0.1039, 0.0012),
    (50.0, 0.1016, 0.0008),
    (70.0, 0.0995, 0.0012),
    (100.0, 0.0990, 0.0012),
];

fn lfit() -> LFit {
    let data: Vec<(f64, f64, f64)> = DATA.iter().map(|(x, y, ey)| (*x, y * x, ey * x)).collect();

    let r = fit_l(&data);
    println!("{r}");
    r
}

// LINEAR FIT: FIT OF INPUT DATA BY STRAIGHT LINE Y=A*X+B.
//
// INPUT:   GAUSSIAN DATA          Y (I), I=1,...,N
//          AND THEIR ERROR BARS   YE(I), I=1,...,N.
//
// OUTPUT:  CONSTANTS  A(1), A(2) WITH ERRORS BARS SGA(1), SGA(2).
//          CORRESPONDING TO THE FIT Y=A(1)+A(2)*X. IF N.GE.3:
//          ALSO  CHI2  AND  THE GOODNESS OF FIT Q ARE RETURNED.

#[derive(Debug, PartialEq)]
struct LFit {
    a: (f64, f64),
    sga: (f64, f64),
    chi2: f64,
    q: f64,
    cov: [[f64; 2]; 2],
}

impl fmt::Display for LFit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "LFit:");
        writeln!(f, "A(1) = {} +/- {}", self.a.0, self.sga.0);
        writeln!(f, "A(2) = {} +/- {}", self.a.1, self.sga.1);
        writeln!(f, "\nCOVARIANCE MATRIX:");
        writeln!(f, "{} {}", self.cov[0][0], self.cov[0][1]);
        writeln!(f, "{} {}", self.cov[1][0], self.cov[1][1]);
        writeln!(f, "\nSTATISTICAL ANALYSIS:");
        writeln!(f, "CHI2 = {}", self.chi2);
        writeln!(f, "   Q = {}", self.q)
    }
}

fn fit_l(data: &[(f64, f64, f64)]) -> LFit {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut s = 0.0;
    for (x, y, sgy) in data {
        let wt = 1.0 / sgy.powi(2);
        s += wt;
        sx += x * wt;
        sy += y * wt;
    }
    let sxos = sx / s;

    let mut stt = 0.0;
    let mut a = (0.0, 0.0);
    for (x, y, sgy) in data {
        let t = (x - sxos) / sgy;
        stt += t * t;
        a.1 += t * y / sgy;
    }
    a.1 /= stt;
    a.0 = (sy - sx * a.1) / s;

    let chi2 = data
        .iter()
        .map(|(x, y, sgy)| ((y - a.0 - a.1 * x) / sgy).powi(2))
        .sum();
    let q = if data.len() > 2 {
        1.0 - gamma_p(0.5 * (data.len() - 2) as f64, 0.5 * chi2)
    } else {
        0.0
    };
    let sga = (
        ((1.0 + sx.powi(2) / (s * stt)) / s).sqrt(),
        (1.0 / stt).sqrt(),
    );
    let cov = -sx / (s * stt);

    LFit {
        a,
        sga,
        chi2,
        q,
        cov: [[sga.0.powi(2), cov], [cov, sga.1.powi(2)]],
    }
}

#[cfg(test)]
mod tests {
    use crate::{lfit, LFit};
    #[test]
    fn lfit_test() {
        let expected = LFit {
            a: (0.17829769783513708, 0.09782177231480599),
            sga: (0.0192082135777534, 0.0007854559197949675),
            chi2: 2.07554626333672,
            q: 0.7218661665309511,
            cov: [
                [0.00036895546884859, -0.000013293884969307767],
                [-0.000013293884969307767, 0.0000006169410019409585],
            ],
        };

        let r = lfit();
        assert_eq!(r, expected);
    }
}
