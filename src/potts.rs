use crate::marsaglia::Marsaglia;

pub struct Potts<const ND: usize, const MS: usize, const NQ: usize> {
    pub ns: usize, // const
    ipf: [[usize; ND]; MS],
    ipb: [[usize; ND]; MS],
    ista: [usize; MS],
    idel: [[usize; NQ]; NQ],
    pub nlink: usize,
    ix: [usize; ND],
    nla: [usize; ND],
    //wrat: Vec<[f64;4*ND+1]>,
}

impl<const ND: usize, const MS: usize, const NQ: usize> Potts<ND, MS, NQ> {
    pub fn new(nla: [usize; ND]) -> Potts<ND, MS, NQ> {
        let ns = nla.iter().fold(1, |acc, x| acc * x); //#lattice sites from edges

        let mut potts = Potts {
            ns,
            ipf: [[0; ND]; MS],
            ipb: [[0; ND]; MS],
            ista: [0; MS],
            idel: [[0; NQ]; NQ],
            nlink: ND * ns,
            ix: [0; ND],
            nla,
        };

        for is in 0..ns {
            for id in 0..ND {
                potts.ixcor(is + 1);
                //Forward (backward) step with periodic bounday conditions:
                potts.ix[id] = (potts.ix[id] + 1) % potts.nla[id];
                potts.ipf[is][id] = potts.calc_is();
            }
            for id in 0..ND {
                potts.ixcor(is + 1);
                //Backward pointer (notice periodic boundary conditions):
                potts.ix[id] = (potts.ix[id] + potts.nla[id] - 1) % potts.nla[id];
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
        for (vix, vnla) in self.ix.iter().zip(self.nla.iter()) {
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
            nspart /= self.nla[id];
            self.ix[id] = (js - 1) / nspart;
        }
    }

    pub fn action(&self) -> usize {
        let mut iact = 0;
        for is in 0..self.ns {
            let ista1 = self.ista[is];
            for id in 0..ND {
                let ista2 = self.ista[self.ipf[is][id] - 1];
                iact += self.idel[ista1][ista2];
            }
        }
        iact
    }

    // Assigns random (i.e. beta=0) values 0,..,nq-1 to the states ista(is).
    pub fn ran(&mut self, rng: &mut Marsaglia) {
        for is in 0..self.ns {
            self.ista[is] = (NQ as f64 * rng.uni()) as usize;
        }
    }

    /// Assigns the value iq to the states ista(is) kept in potts
    /// and returns the corresponding action variable iact.
    pub fn order(&mut self, iq: usize) -> usize {
        for e in self.ista.iter_mut() {
            *e = iq;
        }
        ND * self.ns
    }

    // Initialize ratios of weights for MUCA or canonical Metropolis
    //fn wght_init(&mut self) {
    //    let mut a = [0.0f64;4*ND+1];
    //    for (id,e) in a.iter_mut().enumerate() {
    //        let i = id as f64 - 2.0;
    //        a[id] = (2.0*beta*i).exp();
    //    }
    //    self.wrat = vec![a;nlink];
    //}

    //fn get_wrat(&self, id: i32,ilink: usize) {
    //    self.wrat[ilink][2*ND+id]
    //}
}
