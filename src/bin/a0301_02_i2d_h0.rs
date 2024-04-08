// 2d Ising model: Energy histogram from naive sampling.
use marsaglia_rs::marsaglia::Marsaglia;
use gnuplot::{Caption, Color, Figure, PointSymbol};

const NDAT: usize = 100_000;
const ND: usize = 2;
const ML: usize = 100;
const MS: usize = ML.pow(ND as u32);
//const MLINK: usize = ND * MS;

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
        }
    }

    let mut fg = Figure::new();
    fg.set_title("df all");
    fg.axes2d()
        .lines(&x, &y, &[Caption("Random Sampling"), Color("red")])
		.y_error_bars(&x,&y,&ey, &[Caption(r"y\_error\_bars"), PointSymbol('T'), Color("blue")],);
	

    fg.show().unwrap();
}
