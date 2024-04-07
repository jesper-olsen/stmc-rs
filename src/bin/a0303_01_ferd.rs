//Ising Model on finite lattice after
//FERDINAND and FISHER, PR 185 (1969) 832.

use core::f64::consts::PI;

const NL: usize = 80;
const ML: usize = 80;
const NLL: usize = 2 * NL - 1;
const NBETA: usize = 1;
const DBETA: f64 = 0.4;

fn main() {
    let mut gm = [0.0f64; NLL + 1];
    let mut gmpr = [0.0f64; NLL + 1];
    let mut gm2pr = [0.0f64; NLL + 1];

    for ibeta in 1..=NBETA {
        let beta = DBETA * ibeta as f64;
        gff(beta, &mut gm, &mut gmpr, &mut gm2pr);
        let (_, ueng, _, _, _) = partfn(beta, &gm, &gmpr, &gm2pr);

        println!("2dI_ts_ferdi.d");
        println!("{:10.4}{:10.4}", 0.0, ueng);
        println!("{:10.4}{:10.4}", 200.0, ueng);
    }
}

fn gff(beta: f64, gm: &mut [f64], gmpr: &mut [f64], gm2pr: &mut [f64]) {
    let mut c = [0.0f64; NLL];
    let bk = 2.0 * beta;
    gm[0] = bk + beta.tanh().ln();
    gmpr[0] = 2.0 + 1.0 / beta.tanh() - beta.tanh();
    let sqtnh = beta.tanh().powi(2);
    gm2pr[0] = sqtnh - 1.0 / sqtnh;
    for ill in 1..=NLL {
        c[ill - 1] = bk.cosh() * (1.0 / bk.tanh()) - (PI * ill as f64 / (NL as f64)).cos();
        let csq = (c[ill - 1] * c[ill - 1] - 1.0).sqrt();
        gm[ill] = (c[ill - 1] + csq).ln();

        let taprn = 1.0 - 1.0 / bk.tanh().powi(2);
        let cpr = 2.0 * bk.sinh() / bk.tanh() + 2.0 * bk.cosh() * taprn;
        gmpr[ill] = cpr / csq;

        let c2pr1 = 2.0 * taprn * (bk.sinh() - bk.cosh() / bk.tanh());
        let c2pr = 4.0 * (c2pr1 + bk.cosh() / bk / bk.tanh());
        gm2pr[ill] = (c2pr - gmpr[ill] * gmpr[ill] * c[ill - 1]) / csq;
    }
}

fn partfn(beta: f64, gm: &[f64], gmpr: &[f64], gm2pr: &[f64]) -> (f64, f64, f64, f64, f64) {
    let bk = 2.0 * beta;
    let xml = ML as f64;
    let xnl = NL as f64;

    //  WE PUT IN A MULTIPLICATION FACTOR Z0 TO PREVENT THE PARTITION
    //  FUNCTION BLOW UP.

    let mut zz1 = 0.0;
    let mut zzp1 = 0.0;
    let mut zzpp1 = 0.0;
    let mut zz2 = 0.0;
    let mut zzp2 = 0.0;
    let mut zzpp2 = 0.0;

    let mut zz3 = 0.0;
    let mut zzp3 = 0.0;
    let mut zzpp3 = 0.0;

    let mut zz4 = 0.0;
    let mut zzp4 = 0.0;
    let mut zzpp4 = 0.0;

    let mut z4sign = 1.0;

    for i in 0..NL {
        let ipi = i + i;
        let ipip1 = i + i + 1;
        zz1 += gm[ipip1] * xml / 2.0 + (1.0 + (-gm[ipip1] * xml).exp()).ln();

        zzp1 += gmpr[ipip1] * xml / 2.0 * (1.0 - (-gm[ipip1] * xml).exp())
            / (1.0 + (-gm[ipip1] * xml).exp());
        zzpp1 += gm2pr[ipip1] * xml / 2.0 * (1.0 - (-gm[ipip1] * xml).exp())
            / (1.0 + (-gm[ipip1] * xml).exp())
            + (gmpr[ipip1] * xml / (1.0 + (-gm[ipip1] * xml).exp())).powi(2)
                * (-gm[ipip1] * xml).exp();

        zz2 += gm[ipip1] * xml / 2.0 + (1.0 - (-gm[ipip1] * xml).exp()).ln();
        zzp2 += gmpr[ipip1] * xml / 2.0 * (1.0 + (-gm[ipip1] * xml).exp())
            / (1.0 - (-gm[ipip1] * xml).exp());

        zzpp2 += gm2pr[ipip1] * xml / 2.0 * (1.0 + (-gm[ipip1] * xml).exp())
            / (1.0 - (-gm[ipip1] * xml).exp())
            - (gmpr[ipip1] * ML as f64 / (1.0 - (-gm[ipip1] * xml).exp())).powi(2)
                * (-gm[ipip1] * xml).exp();
        if gm[ipi] >= 0.0 {
            zz3 += gm[ipi] * xml / 2.0 + (1.0 + (-gm[ipi] * xml).exp()).ln();
            zzp3 += gmpr[ipi] * xml / 2.0 * (1.0 - (-gm[ipi] * xml).exp())
                / (1.0 + (-gm[ipi] * xml).exp());
            zzpp3 += gm2pr[ipi] * xml / 2.0 * (1.0 - (-gm[ipi] * xml).exp())
                / (1.0 + (-gm[ipi] * xml).exp())
                + (gmpr[ipi] * xml / (1.0 + (-gm[ipi] * xml).exp())).powi(2)
                    * (-gm[ipi] * xml).exp();

            zz4 += gm[ipi] * ML as f64 / 2.0 + (1.0 - (-gm[ipi] * xml).exp()).ln();
            zzp4 += gmpr[ipi] * xml / 2.0 * (1.0 + (-gm[ipi] * xml).exp())
                / (1.0 - (-gm[ipi] * xml).exp());
            zzpp4 += gm2pr[ipi] * xml / 2.0 * (1.0 + (-gm[ipi] * xml).exp())
                / (1.0 - (-gm[ipi] * xml).exp())
                - (gmpr[ipi] * xml / (1.0 - (-gm[ipi] * xml).exp())).powi(2)
                    * (-gm[ipi] * xml).exp();
        } else {
            zz3 += -gm[ipi] * xml / 2.0 + (1.0 + (gm[ipi] * xml).exp()).ln();
            zzp3 += -gmpr[ipi] * xml / 2.0 * (1.0 - (gm[ipi] * xml).exp())
                / (1.0 + (gm[ipi] * xml).exp());
            zzpp3 += -gm2pr[ipi] * xml / 2.0 * (1.0 - (gm[ipi] * xml).exp())
                / (1.0 + (gm[ipi] * xml).exp())
                + (gmpr[ipi] * xml / (1.0 + (gm[ipi] * xml).exp())).powi(2) * (gm[ipi] * xml).exp();
            zz4 += -gm[ipi] * xml / 2.0 + (1.0 - (gm[ipi] * xml).exp()).ln();
            z4sign = -z4sign;
            zzp4 += -gmpr[ipi] * xml / 2.0 * (1.0 + (gm[ipi] * xml).exp())
                / (1.0 - (gm[ipi] * xml).exp());
            zzpp4 += -gm2pr[ipi] * xml / 2.0 * (1.0 + (gm[ipi] * xml).exp())
                / (1.0 - (gm[ipi] * xml).exp())
                - (gmpr[ipi] * xml / (1.0 - (gm[ipi] * xml).exp())).powi(2) * (gm[ipi] * xml).exp();
        }
    }

    let mut zl =
        zz1 + (1.0 + (zz2 - zz1).exp() + (zz3 - zz1).exp() + z4sign * (zz4 - zz1).exp()).ln();
    let zlp = zzp1 * (zz1 - zl).exp()
        + zzp2 * (zz2 - zl).exp()
        + zzp3 * (zz3 - zl).exp()
        + z4sign * zzp4 * (zz4 - zl).exp();
    let ueng = -1.0 / bk.tanh() - zlp / xml / xnl;
    zzpp1 += zzp1.powi(2);
    zzpp2 += zzp2.powi(2);
    zzpp3 += zzp3.powi(2);
    zzpp4 += zzp4.powi(2);
    let cht1 = zzpp1 * (zz1 - zl).exp()
        + zzpp2 * (zz2 - zl).exp()
        + zzpp3 * (zz3 - zl).exp()
        + z4sign * zzpp4 * (zz4 - zl).exp();
    let cht = (cht1 - zlp.powi(2)) * (bk * bk) / (4.0 * xnl * xml) - (bk / bk.sinh()).powi(2) / 2.0;
    zl += xml * xnl / 2.0 * (2.0 * bk.sinh()).ln() - (2.0f64).ln();

    let f = -zl / (xml * xnl) / beta;
    let s = (ueng - f) * beta;
    (f, ueng, s, cht, zl)
}

//fn sign(x: f64, y: f64) -> f64 {
//    if y > 0.0 {
//        x
//    } else if y < 0.0 {
//        -x
//    } else {
//        0.0
//    }
//}
