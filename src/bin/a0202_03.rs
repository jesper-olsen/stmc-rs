
use marsaglia_rs::marsaglia::Marsaglia;
use marsaglia_rs::steb::{steb0};
use marsaglia_rs::{qtiles};

/// Variance fractiles: experimental verification.
fn main() {
    const NDAT:usize=10000;
    const ND2:usize=2;
    let mut var = [0.0f64; NDAT];
    let mut data = [0.0f64; NDAT];
    let mut rng = Marsaglia::new(12, 34, 56, 78);
    for i in 0..NDAT {
        data[0]=rng.gauss();
        data[1]=rng.gauss();
        (_,var[i])=steb0(&data[0..ND2]);
        var[i]=2.0*var[i].powi(2);
    }
    var.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q=0.15;
    let (vq1,vq2) = qtiles(&var,q).unwrap(); 
    println!("q: {q} vq1: {vq1} vq2: {vq2}");
}

//      PARAMETER (IUO=6,IUD=10,NDAT=10000,ND2=2,ISEED1=1,ISEED2=0)
//      DIMENSION VAR(NDAT),DATA(ND2)
//C
//      CALL RMASET(IUO,IUD,ISEED1,ISEED2,'nofile')
//      WRITE(IUO,*) ' '
//      DO N=1,NDAT
//        CALL RMAGAU(DATA(1),DATA(2))
//        CALL STEB0(ND2,DATA,DM,DV,VAR(N))
//        VAR(N)=2.0D00*VAR(N)**2
//      END DO
//C
//      CALL HEAPSORT(NDAT,VAR)
//      Q=0.15D00
//      CALL QTILES(NDAT,VAR,Q,VQ1,VQ2)
//      WRITE(IUO,100) Q,VQ1,VQ2
//100   FORMAT(1X,'Q,VQ1,VQ2: ',3F18.8)
//C
