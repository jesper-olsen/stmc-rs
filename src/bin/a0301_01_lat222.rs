// Bookkeeping for a hypercubic 2x2x2 lattice (lat

const ND: usize = 3;
const ML: usize = 3;
const MS: usize = ML.pow(ND as u32);
const NLA: [usize; ND] = [2; ND];

struct Lat {
    ns: usize,
    ipf: [[usize; ND]; MS],
    ipb: [[usize; ND]; MS],
    ix: [usize; ND],
}

impl Lat {
    pub fn new() -> Lat {
        let ns = NLA.iter().fold(1, |acc, x| acc * x);
        let mut lat = Lat {
            ns,
            ipf: [[0; ND]; MS],
            ipb: [[0; ND]; MS],
            ix: [0; ND],
        };
        for is in 1..=ns {
            for id in 0..ND {
                lat.ixcor(is);
                //Forward (backward) step with periodic bounday conditions:
                lat.ix[id] = (lat.ix[id] + 1) % NLA[id];
                lat.ipf[is][id] = lat.calc_is();
            }
            for id in 0..ND {
                lat.ixcor(is);
                //Backward pointer (notice periodic boundary conditions):
                lat.ix[id] = (lat.ix[id] + NLA[id] - 1) % NLA[id];
                lat.ipb[is][id] = lat.calc_is();
            }
        }
        lat
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

    fn calc_is(&self) -> usize {
        let mut is = 1;
        let mut nsa = 1;
        for (vix, vnla) in self.ix.iter().zip(NLA.iter()) {
            is += vix * nsa;
            nsa *= *vnla;
        }
        is
    }
}

fn main() {
    println!("   is   ix(1) (2) (3)   ipf(is,1) (,2) (,3)   ipb(is,1) (,2) (,3)");
    let mut lat = Lat::new();

    for is in 1..=lat.ns {
        lat.ixcor(is);
        println!(
            " {is:4}   {:?}       {:?}             {:?}",
            lat.ix, &lat.ipf[is], &lat.ipb[is]
        );
    }
}
