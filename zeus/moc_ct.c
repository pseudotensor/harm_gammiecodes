
/* step z magnetic field, velocity using MOC, CT */

#include "decs.h"

void moc_ct()
{
	double **bxstar,**bystar,**vxstar,**vystar ;
	double **dqv,**dqb,**bxprim,**byprim ;
	double **emf,bxa,bya ; 
	double Dx1,Dx2,rhoa,vaa,vxa,vya,srhoa,vp,vm,bp,bm,Dy1,Dy2 ;
	double sgn_va ;
	double pr ;
	int i,j ;
	void bound_var(double **var) ;
	void dqx_calc(double **var, double **dq) ;
	void dqy_calc(double **var, double **dq) ;

	dqv = work3 ;
	dqb = work4 ;
	vxstar = work5 ;
	vystar = work6 ;
	bxstar = work7 ;
	bystar = work8 ;
	bxprim = work9 ;
	byprim = work10 ;

	/* sweep in x-direction */
	/* first get slopes */

	dqx_calc(vy,dqv) ;
	dqx_calc(By,dqb) ;

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		rhoa = 0.25*(rho[i][j] + rho[i-1][j] +
			rho[i][j-1] + rho[i-1][j-1]) ;
		srhoa = sqrt(rhoa) ;
		vaa = 0.5*((Bx[i][j]+Bx[i][j-1])/srhoa) ;
		sgn_va = copysign(1.,vaa) ;
		vaa = fabs(vaa) ;
		vxa = 0.5*(vx[i][j]+vx[i][j-1]) ;
		Dx1 = dt*vaa/(a*dx) ;
		Dx2 = dt*vxa/(a*dx) ;

		if(Dx1 + Dx2 > 0.) {
			vp = vy[i-1][j] + (1. - Dx1 - Dx2)*dqv[i-1][j] ;
			bp = By[i-1][j] + (1. - Dx1 - Dx2)*dqb[i-1][j] ;
		}
		else {
			vp = vy[i][j] + (-1. - Dx1 - Dx2)*dqv[i][j] ;
			bp = By[i][j] + (-1. - Dx1 - Dx2)*dqb[i][j] ;
		}
		if(-Dx1 + Dx2 > 0.) {
			vm = vy[i-1][j] + (1. + Dx1 - Dx2)*dqv[i-1][j] ;
			bm = By[i-1][j] + (1. + Dx1 - Dx2)*dqb[i-1][j] ;
		}
		else {
			vm = vy[i][j] + (-1. + Dx1 - Dx2)*dqv[i][j] ;
			bm = By[i][j] + (-1. + Dx1 - Dx2)*dqb[i][j] ;
		}
		
                /* solution to Stone & Norman eqtn 43,44-- use constant rho */
                vystar[i][j] = 0.5*(vm + vp + sgn_va*(bm - bp)/srhoa) ;
                bystar[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*srhoa) ;

		/* find bprim, for updating velocities */
                vp = vy[i-1][j] + (1. - Dx1)*dqv[i-1][j] ;
                bp = By[i-1][j] + (1. - Dx1)*dqb[i-1][j] ;
                vm = vy[i][j] + (-1. + Dx1)*dqv[i][j] ;
                bm = By[i][j] + (-1. + Dx1)*dqb[i][j] ;

                byprim[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*srhoa) ;
	}
	bound_var(vystar) ;
	bound_var(bystar) ;
	bound_var(byprim) ;

	/* sweep in y-direction */
	/* first get slopes */
	dqy_calc(vx,dqv) ;
	dqy_calc(Bx,dqb) ;

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		rhoa = 0.25*(rho[i][j] + rho[i-1][j] +
			rho[i][j-1] + rho[i-1][j-1]) ;
		srhoa = sqrt(rhoa) ;
		vaa = 0.5*((By[i][j]+By[i-1][j])/srhoa) ;
		sgn_va = copysign(1.,vaa) ;
		vaa = fabs(vaa) ;
		vya = 0.5*(vy[i][j]+vy[i-1][j]) ;
		Dy1 = dt*vaa/(a*dy) ;
		Dy2 = dt*vya/(a*dy) ;

                /* values at the foot of the plus characteristic
                  are, by convention, in i-1 zone */
		if(Dy1 + Dy2 > 0.) {
			vp = vx[i][j-1] + (1. - Dy1 - Dy2)*dqv[i][j-1] ;
			bp = Bx[i][j-1] + (1. - Dy1 - Dy2)*dqb[i][j-1] ;
		}
		else {
			vp = vx[i][j] + (-1. - Dy1 - Dy2)*dqv[i][j] ;
			bp = Bx[i][j] + (-1. - Dy1 - Dy2)*dqb[i][j] ;
		}
		if(-Dy1 + Dy2 > 0.) {
			vm = vx[i][j-1] + (1. + Dy1 - Dy2)*dqv[i][j-1] ;
			bm = Bx[i][j-1] + (1. + Dy1 - Dy2)*dqb[i][j-1] ;
		}
		else {
			vm = vx[i][j] + (-1. + Dy1 - Dy2)*dqv[i][j] ;
			bm = Bx[i][j] + (-1. + Dy1 - Dy2)*dqb[i][j] ;
		}
		
                /* solution to Stone & Norman eqtn 43,44-- use constant rho */
                vxstar[i][j] = 0.5*(vm + vp + sgn_va*(bm - bp)/srhoa) ;
                bxstar[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*srhoa) ;

		/* find bprim, to update velocities */
                vp = vx[i][j-1] + (1. - Dy1)*dqv[i][j-1] ;
                vm = vx[i][j] + (-1. + Dy1)*dqv[i][j] ;
                bp = Bx[i][j-1] + (1. - Dy1)*dqb[i][j-1] ;
                bm = Bx[i][j] + (-1. + Dy1)*dqb[i][j] ;

                bxprim[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*srhoa) ;
	}
	bound_var(vxstar) ;
	bound_var(bxstar) ;
	bound_var(bxprim) ;

	/* calculate emf */
	emf = work1 ;
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		emf[i][j] = vxstar[i][j]*bystar[i][j] - vystar[i][j]*bxstar[i][j] ;
	}
	bound_var(emf) ;

	/* update field */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		Bx[i][j] += dt*(emf[i][j+1] - emf[i][j])/(a*dy) ;
		By[i][j] += dt*(emf[i][j] - emf[i+1][j])/(a*dx) ;
	}
	bound_var(Bx) ;
	bound_var(By) ;

	/* update velocity */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		rhoa = 0.5*(rho[i][j]+rho[i-1][j]) ;
		bya = 0.25*(By[i][j] + By[i-1][j] + 
			By[i][j+1] + By[i-1][j+1]) ;

		vx[i][j] += dt*bya*(bxprim[i][j+1] - bxprim[i][j])/(a*dy*rhoa) ;

		rhoa = 0.5*(rho[i][j]+rho[i][j-1]) ;
		bxa = 0.25*(Bx[i][j] + Bx[i+1][j] + 
			Bx[i][j-1] + Bx[i+1][j-1]) ;

		vy[i][j] += dt*bxa*(byprim[i+1][j] - byprim[i][j])/(a*dx*rhoa) ;
	}
	bound_var(vx) ;
	bound_var(vy) ;

	/* done! */
}
