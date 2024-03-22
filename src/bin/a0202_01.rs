use marsaglia_rs::chi2::{chi2_df, chi2_xq};

const NMAX: usize = 8;

fn main() {
    let mut f = [0.0f64; NMAX];
    let q = [0.025, 0.15, 0.50, 0.85, 0.975];
    let mut xq = [0.0f64; 5];
    let mut nf: usize;

    println!("\n CHI2 distribution function\n");
    print!("   CHI2 \\ N ");
    for n in 1..=NMAX {
        print!(" {n:4}  ");
    }
    println!();
    for i in 1..=NMAX {
        let chi2 = i as f64;
        print!("   {chi2:1.2}    ");
        for n in 0..NMAX {
            f[n] = chi2_df(chi2, n + 1);
            print!(" {:2.3} ", f[n]);
        }
        println!();
    }

    println!("\n Fractiles for the CHI2 distribution\n");
    print!("   N \\ q ");
    for n in 0..q.len() {
        print!(" {: >8.3}  ", q[n]);
    }
    println!("\n");
    for n in 1..=NMAX {
        let nf = n;
        print!("   {n}     ");
        for i in 0..q.len() {
            xq[i] = chi2_xq(q[i], n);
            print!(" {: >8.3}  ", xq[i]);
        }
        println!();
    }

    println!("\n CHI2 per degree of freedom distribution function");
    print!("   CHI2 \\ N ");
    for n in 1..=NMAX {
        print!(" {n:4}  ");
    }

    println!("\n Fractiles for the CHI2 per degree of freedom distribution\n");
    print!("   N \\ q ");
    for n in 0..q.len() {
        print!(" {: >8.3}  ", q[n]);
    }
}

//      WRITE(IUO,*) ' '
//      WRITE(IUO,*)
//     & ' CHI2 per degree of freedom distribution function '
//      WRITE(IUO,*) ' '
//      WRITE(IUO,102) (N,N=1,NMAX)
//102   FORMAT(1X,'(CHI2/N) \\ N',3X,10(1X,I4,2X)) ! UNIX
//C 102   FORMAT(1X,'(CHI2/N) \ N',3X,10(1X,I4,2X)) ! DOS
//      WRITE(IUO,*) ' '
//      DO I=1,NMAX
//C CHI2 PER DEGREE OF FREEDOM:
//        X=(I*ONE)/FOUR
//C LOOP OVER DEGREES OF FREEDOM:
//        DO N=1,NMAX
//          NF=N
//          F(N)=CHI2PDF_DF(X)
//        END DO
//        WRITE(IUO,100) X,F
//      END DO
//      WRITE(IUO,*) ' '
//C
//      WRITE(IUO,*) ' '
//      WRITE(IUO,*)
//     & ' Fractiles for the CHI2 per degree of freedom distribution'
//      WRITE(IUO,*) ' '
//      WRITE(IUO,103) (Q(I),I=1,NQ)
//      WRITE(IUO,*) ' '
//C LOOP OVER DEGREES OF FREEDOM:
//      DO N=1,NMAX
//        NF=N
//C LOOP OVER Q-TILES:
//        DO I=1,NQ
//          XQ(I)=CHI2PDF_XQ(Q(I))
//        END DO
//        WRITE(IUO,200) N,XQ
//      END DO
//      WRITE(IUO,*) ' '
//C
//      STOP
//      END
//
//      INCLUDE '../../ForLib/chi2_df.f'
//      INCLUDE '../../ForLib/chi2_xq.f'
//      INCLUDE '../../ForLib/chi2pdf_df.f'
//      INCLUDE '../../ForLib/chi2pdf_xq.f'
//      INCLUDE '../../ForLib/fi1.f'
//      INCLUDE '../../ForLib/gamma_p.f'
//      INCLUDE '../../ForLib/gamma_ln.f'
