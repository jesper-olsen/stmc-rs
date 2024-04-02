// Bookkeeping for a hypercubic 2x2x2 lattice (lat

const ND: usize = 3;
const ML: usize = 3;
const MS: usize = ML.pow(ND as u32);
const NLA: [usize; ND] = [2; ND];
struct Lat {
    ns: usize,
    ipf: [[i32; ND]; MS],
    ipb: [[i32; ND]; MS],
    ix: [i32; ND],
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
                lat.ixcor(is as i32);
                //Forward (backward) step with periodic bounday conditions:
                lat.ix[id] = (lat.ix[id] + 1) % NLA[id] as i32;
                lat.ipf[is][id] = lat.calc_is();
            }
            for id in 0..ND {
                lat.ixcor(is as i32);
                //Backward pointer (notice periodic boundary conditions):
                lat.ix[id] = (lat.ix[id] - 1 + NLA[id] as i32) % NLA[id] as i32;
                lat.ipb[is][id] = lat.calc_is();
            }
        }
        lat
    }

    fn ixcor(&mut self, is: i32) {
        let ns = NLA.iter().fold(1, |acc, x| acc * x);
        let mut nspart = ns as i32;
        let mut js = is;
        for id in (0..ND).rev() {
            if id < ND - 1 {
                js -= self.ix[id + 1] * nspart ;
            }
            nspart /= NLA[id] as i32;
            self.ix[id] = (js - 1) / nspart;
        }
    }

    fn calc_is(&self) -> i32 {
        let mut is = 1;
        let mut nsa = 1;
        for (vix, vnla) in self.ix.iter().zip(NLA.iter()) {
            is += vix * nsa;
            nsa *= *vnla as i32;
        }
        is
    }
}

fn main() {
    println!("   is   ix(1) (2) (3)   ipf(is,1) (,2) (,3)   ipb(is,1) (,2) (,3)");
    let mut lat = Lat::new();

    for is in 1..=lat.ns {
        lat.ixcor(is as i32);
        let ipf = &lat.ipf[is];
        let ipb = &lat.ipb[is];
        println!(" {is:4}   {:?}       {:?}             {:?}", lat.ix, ipf, ipb);
    }
}
