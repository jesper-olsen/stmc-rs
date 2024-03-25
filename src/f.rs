use crate::beta::beta_i;
use crate::gamma::gamma_ln;
use crate::chi2::fi1;

/// VARIANCE RATIO DISTRIBUTION FUNCTION.
pub fn f_df(f: f64, nf1: u32, nf2: u32) -> f64 {
    let nf1 = nf1 as f64;
    let nf2 = nf2 as f64;
    let xf1h=0.5*nf1;
    let xf2h=0.5*nf2;
    let x = nf2 / (nf1*f+nf2);
    1.0-beta_i(x,xf2h,xf1h)
}

// VARIANCE RATIO TEST (F-TEST): COMPARISION OF TWO VARIANCES.
pub fn ftest(eb1: f64, eb2: f64, ndat1: u32, ndat2: u32) -> f64 {
    let nf1 = ndat1-1; // degrees of freedom
    let nf2 = ndat2-1;
    let va1 = eb1.powi(2)*ndat1 as f64;
    let va2 = eb2.powi(2)*ndat2 as f64;
    let f = va1/va2; //  Definitions take per degree of freedom already into account!
    let q=2.0*f_df(f, nf1, nf2);
    if q>1.0 {
        2.0-q
    } else {
        q
    }
}

//C VARIANCE RATIO PROBABILITY DENSITY
pub fn f_pd(f: f64, nf1: u32, nf2: u32) -> f64 {
    let f1h = 0.5*nf1 as f64;
    let f2h = 0.5*nf2 as f64;
    let y = f1h*f/f2h;
    let fh = f1h+f2h;
    f1h*-(gamma_ln(f1h)-gamma_ln(f2h)+gamma_ln(fh)).exp()*y.powf(f1h-1.0)*(y+1.0).powf(-fh)/f2h
}

// VARIANCE RATIO DISTRIBUTION FUNCTION.
pub fn f_qdf(f: f64, nf1: u32, nf2: u32) -> f64 {
    let nf1 = nf1 as f64;
    let nf2 = nf2 as f64;
    let xf1h=0.5*nf1;
    let xf2h=0.5*nf2;
    let x = nf2 as f64 / (nf1*f+nf2);
    let f_qdf=1.0-beta_i(x,xf2h,xf1h);
    if f_qdf > 0.5 {
        1.0-f_qdf
    } else {
        f_qdf
    }
}

pub fn f_xq(q: f64, nf1: u32, nf2: u32) -> f64 {
    if q> 0.5 {
        1.0
    } else {
        let x1;
        let mut x2=1.0+1.0/10.0;
        if q>0.5 {
            x1 = 1.0-1.0/10.0;
            loop {
                x2 += 1.0;
                if f_df(x2, nf1, nf2)> q {
                    break;
                }
           }
        } else {
            x1=0.0;
        }
        let f = |x: f64| f_df(x, nf1, nf2);
        fi1(f, q, x1, x2)
    }
} 
