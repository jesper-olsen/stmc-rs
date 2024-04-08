// 2d Ising model: Energy histogram from naive sampling.
use gnuplot::{AxesCommon, Caption, Color, Figure, PointSymbol, AutoOption};
use marsaglia_rs::marsaglia::Marsaglia;
use std::fs::File;
use std::io::Write;

const NDAT: usize = 100_000;
const ND: usize = 2;
const ML: usize = 100;
const MS: usize = ML.pow(ND as u32);
const MLINK: usize = ND * MS;

const NQ: usize = 2;
const NLA: [usize; ND] = [20, 20];

struct Potts {
    ns: usize, // const
    ipf: [[usize; ND]; MS],
    ipb: [[usize; ND]; MS],
    ista: [usize; MS],
    idel: [[usize; NQ]; NQ],
    nlink: usize,
    ix: [usize; ND],
}

impl Potts {
    pub fn new() -> Potts {
        let ns = NLA.iter().fold(1, |acc, x| acc * x); //#lattice sites from edges

        let mut potts = Potts {
            ns,
            ipf: [[0; ND]; MS],
            ipb: [[0; ND]; MS],
            ista: [0; MS],
            idel: [[0; NQ]; NQ],
            nlink: ND * ns,
            ix: [0; ND],
        };

        for is in 0..ns {
            for id in 0..ND {
                potts.ixcor(is + 1);
                //Forward (backward) step with periodic bounday conditions:
                potts.ix[id] = (potts.ix[id] + 1) % NLA[id];
                potts.ipf[is][id] = potts.calc_is();
            }
            for id in 0..ND {
                potts.ixcor(is + 1);
                //Backward pointer (notice periodic boundary conditions):
                potts.ix[id] = (potts.ix[id] + NLA[id] - 1) % NLA[id];
                potts.ipb[is][id] = potts.calc_is();
            }
        }
        for i in 0..NQ {
            potts.idel[i][i] = 1;
        }
        potts
    }

    fn calc_is(&self) -> usize {
        let mut is = 1;
        let mut nsa = 1;
        for (vix, vnla) in self.ix.iter().zip(NLA.iter()) {
            is += vix * nsa;
            nsa *= *vnla;
        }
        is
    }

    fn ixcor(&mut self, is: usize) {
        let mut nspart = self.ns;
        let mut js = is;
        for id in (0..ND).rev() {
            if id < ND - 1 {
                js -= self.ix[id + 1] * nspart;
            }
            nspart /= NLA[id];
            self.ix[id] = (js - 1) / nspart;
        }
    }

    fn action(&self) -> usize {
        let mut iact = 0;
        for is in 0..self.ns {
            // MS?
            let ista1 = self.ista[is];
            for id in 0..ND {
                let ista2 = self.ista[self.ipf[is][id] - 1];
                iact += self.idel[ista1][ista2];
            }
        }
        iact
    }

    // Assigns random (i.e. beta=0) values 0,..,nq-1 to the states ista(is).
    fn ran(&mut self, rng: &mut Marsaglia) {
        for is in 0..self.ns {
            self.ista[is] = (NQ as f64 * rng.uni()) as usize;
        }
    }
}

fn main() {
    let mut rng = Marsaglia::new(12, 34, 56, 78);

    let mut potts = Potts::new();
    let mut ha = vec![0; potts.nlink + 1];
    for idat in 1..=NDAT {
        potts.ran(&mut rng);
        let iact = potts.action();
        ha[iact] += 1;
        let i_e = 2 * iact as i32 - potts.nlink as i32;
        if idat % (NDAT / 10) == 0 {
            println!(" idat,iact,iE: {idat:10}{iact:10}{i_e:10}");
        }
    }

    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut ey = Vec::new();
    let mut v = Vec::new();
    for iact in 0..=potts.nlink {
        let i_e = potts.nlink as i32 - 2 * iact as i32;
        let es = i_e as f64 / potts.ns as f64;
        let p = ha[iact] as f64 / NDAT as f64;
        let pe = (p * (1.0 - p) / (NDAT - 1) as f64).sqrt(); // Binomial error bar.
        let he = NDAT as f64 * pe;
        if ha[iact] > 0 {
            //println!("{i_e:6}{es:12.5}{p:12.5}{pe:12.5}{h:12.5}{he:12.5}", h=ha[iact] as f64);
            x.push(es);
            y.push(ha[iact] as f64);
            ey.push(he);
            v.push((i_e, es, p, pe, ha[iact] as f64, he));
        }
    }

    let mut fg = Figure::new();
    fg.set_title("df all");
    fg.axes2d()
        .set_x_label("e_s", &[])
        .set_y_label("H", &[])
        .lines(&x, &y, &[Caption("Random Sampling"), Color("red")])
        .y_error_bars(
            &x,
            &y,
            &ey,
            &[Caption(r"y\_error\_bars"), PointSymbol('T'), Color("blue")],
        );
    fg.show().unwrap();

    let mut hae = [0.0; MLINK];
    let mut hb = [0.0; MLINK];
    let mut hbe = [0.0; MLINK];

    let mut iamin = 0;
    let mut iact = 0;
    for (idat, (i_e, _es, _p, _pe, h, he)) in v.iter().enumerate() {
        iact = (-(i_e - potts.nlink as i32) / 2) as usize;
        if idat == 0 {
            iamin = iact;
        }
        ha[iact] = *h as usize;
        hae[iact] = *he;
    }
    let iamax = iact;

    // Re-weighting of the histogram ha() and its error bars hae() at beta0
    //              to the histogram hb() and its error bars hbe() at beta.
    let beta = 0.2f64;
    let beta0 = 0.0;
    let mut hasum = 0.0;
    let mut act0m = 0.0;
    for (iact,e) in ha.iter().enumerate() {
        hasum += *e as f64;
        act0m += (iact * *e) as f64;
    }
    act0m /= hasum;

    let mut hbsum = 0.0;
    for iact in 0..=potts.nlink {
        hb[iact] = 0.0;
        hbe[iact] = 0.0;
        if ha[iact] > 0 {
            let a = (beta - beta0) * 2.0 * (iact as f64 - act0m);
            hb[iact] = ((ha[iact] as f64).ln() + a).exp();
            hbe[iact] = (hae[iact].ln() + a).exp();
            hbsum += hb[iact];
        } // else hb[iact] and ha[iact] are 0
    }

    let factor = hasum / hbsum;
    for iact in 0..=potts.nlink {
        hb[iact] *= factor;
        hbe[iact] *= factor;
    }

    println!("I2d_Hb.d");
    let mut file = File::create("I2d_Hb.d").unwrap();
    let mut x2 = Vec::new();
    let mut y2 = Vec::new();
    let mut ey2 = Vec::new();
    for iact in (iamin..=iamax).step_by(2) {
        let i_e = -2 * iact as i32 + potts.nlink as i32;
        let es = i_e as f64 / potts.ns as f64;
        x2.push(es);
        y2.push(hb[iact]);
        ey2.push(hbe[iact]);
        //println!("{i_e:16.5} {es:12.5} {:12.5} {:12.5}", hb[iact], hbe[iact]);
        writeln!(file, "{i_e:6}{es:12.5}{:12.5}{:12.5}", hb[iact], hbe[iact]).unwrap();
    }

    let mut fg = Figure::new();
    fg.set_title("df all");
    fg.axes2d()
        //.set_x_label("X Axis Label", &[LabelOption::Font("Helvetica"), LabelOption::FontSize(12)])
        .set_x_label("e_s", &[])
        .set_x_range(AutoOption::Fix(0.4), AutoOption::Fix(-2.0))
        .set_y_range(AutoOption::Fix(0.0), AutoOption::Fix(6000.0))
        .set_y_label("H", &[])
                .lines(&x, &y, &[Caption("Random Sampling"), Color("red")])
                .y_error_bars(
                    &x,
                    &y,
                    &ey,
                    &[Caption(r"y\_error\_bars"), PointSymbol('T'), Color("blue")],
                )
        .lines(
            &x2,
            &y2,
            &[Caption("Weighted to \u{03B2}=0.2"), Color("green")],
        )
        .y_error_bars(
            &x2,
            &y2,
            &ey2,
            &[
                Caption(r"y\_error\_bars"),
                PointSymbol('T'),
                Color("orange"),
            ],
        );
    fg.show().unwrap();
}
