
/* step z magnetic field, velocity using MOC */

#include "decs.h"

void step_bz()
{
	void bzsweepx() ;
	void bzsweepy() ;
	static int nsteps = 0 ;
	void spec_diag() ;

	if(nsteps%2 == 0) {
		bzsweepx() ;
		bzsweepy() ;
	}
	else {
		bzsweepy() ;
		bzsweepx() ;
	}

	nsteps++ ;

}

void bzsweepx()
{

	static double **bstar,**vstar,**dqv,**dqb ;
	double Dx,bm,bp,vm,vp ;
	double rhoa,sgn_va ;
	double dqvm,dqvp,dqbm,dqbp,pr,va ;
	int i,j ;
	void bound_var() ;

	vstar = work1 ;	
	bstar = work2 ;
	dqv = work3 ;
	dqb = work4 ;

	/* find vanleer slopes for vz,bz */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		dqvm = (vz[i][j] - vz[i-1][j])/(a*dx) ;
		dqvp = (vz[i+1][j] - vz[i][j])/(a*dx) ;
		dqbm = (Bz[i][j] - Bz[i-1][j])/(a*dx) ;
		dqbp = (Bz[i+1][j] - Bz[i][j])/(a*dx) ; 

		pr = dqvm*dqvp ;
		if(pr > 0.) dqv[i][j] = pr/(dqvm+dqvp) ;
		else dqv[i][j] = 0. ;

		pr = dqbm*dqbp ;
		if(pr > 0.) dqb[i][j] = pr/(dqbm+dqbp) ;
		else dqb[i][j] = 0. ;
	}
	bound_var(dqv) ;
	bound_var(dqb) ;

	/* vstar,bstar, located at zone boundary */
	/* bz,vz: first calculate vstar, bstar */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		rhoa = 0.5*(rho[i][j]+rho[i-1][j]) ;
		va = Bx[i][j]/sqrt(rhoa) ;
		sgn_va = copysign(1.,va) ;
		va = fabs(va) ;
		Dx = va*dt ;

		/* values at the foot of the plus characteristic
		  are, by convention, in i-1 zone */
		vp = vz[i-1][j] + (a*dx - Dx)*dqv[i-1][j] ;
		vm = vz[i][j] - (a*dx - Dx)*dqv[i][j] ;
		bp = Bz[i-1][j] + (a*dx - Dx)*dqb[i-1][j] ;
		bm = Bz[i][j] - (a*dx - Dx)*dqb[i][j] ;

		rhoa = sqrt(rhoa) ;

		/* solution to Stone & Norman eqtn 43,44-- use constant rho */
		vstar[i][j] = 0.5*(vm + vp + sgn_va*(bm - bp)/rhoa) ;
		bstar[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*rhoa) ;
	}
	bound_var(vstar) ;
	bound_var(bstar) ;

	/* step bz,vz */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vz[i][j] += dt*0.5*(Bx[i][j]+Bx[i+1][j])*
			(bstar[i+1][j] - bstar[i][j])/(rho[i][j]*a*dx) ;
		Bz[i][j] += dt*0.5*(Bx[i][j]+Bx[i+1][j])*
			(vstar[i+1][j] - vstar[i][j])/(a*dx) ;
	}
	bound_var(vz) ;
	bound_var(Bz) ;
}

void bzsweepy()
{

	static double **bstar,**vstar,**dqv,**dqb ;
	double Dy,bm,bp,vm,vp ;
	double rhoa,sgn_va ;
	double dqvm,dqvp,dqbm,dqbp,pr,va ;
	int i,j ;
	void bound_var() ;

	vstar = work1 ;	
	bstar = work2 ;
	dqv = work3 ;
	dqb = work4 ;

	/* find vanleer slopes for vz,bz */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		dqvm = (vz[i][j] - vz[i][j-1])/(a*dy) ;
		dqvp = (vz[i][j+1] - vz[i][j])/(a*dy) ;
		dqbm = (Bz[i][j] - Bz[i][j-1])/(a*dy) ;
		dqbp = (Bz[i][j+1] - Bz[i][j])/(a*dy) ;

		pr = dqvm*dqvp ;
		if(pr > 0.) dqv[i][j] = pr/(dqvm+dqvp) ;
		else dqv[i][j] = 0. ;

		pr = dqbm*dqbp ;
		if(pr > 0.) dqb[i][j] = pr/(dqbm+dqbp) ;
		else dqb[i][j] = 0. ;
	}
	bound_var(dqv) ;
	bound_var(dqb) ;

	/* vstar,bstar, located at zone boundary */
	/* bz,vz: first calculate vstar, bstar */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		rhoa = 0.5*(rho[i][j]+rho[i][j-1]) ;
		va = By[i][j]/sqrt(rhoa) ;
		sgn_va = copysign(1.,va) ;
		va = fabs(va) ;
		Dy = va*dt ;

		/* values at the foot of the plus characteristic
		  are, by convention, in j-1 zone */
		vp = vz[i][j-1] + (a*dy - Dy)*dqv[i][j-1] ;
		vm = vz[i][j] - (a*dy - Dy)*dqv[i][j] ;
		bp = Bz[i][j-1] + (a*dy - Dy)*dqb[i][j-1] ;
		bm = Bz[i][j] - (a*dy - Dy)*dqb[i][j] ;

		rhoa = sqrt(rhoa) ;

		/* solution to Stone & Norman eqtn 43,44-- use constant rho */
		vstar[i][j] = 0.5*(vm + vp + sgn_va*(bm - bp)/rhoa) ;
		bstar[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*rhoa) ;
	}
	bound_var(vstar) ;
	bound_var(bstar) ;

	/* step bz,vz */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vz[i][j] += dt*0.5*(By[i][j]+By[i][j+1])*
			(bstar[i][j+1] - bstar[i][j])/(rho[i][j]*a*dy) ;
		Bz[i][j] += dt*0.5*(By[i][j]+By[i][j+1])*
			(vstar[i][j+1] - vstar[i][j])/(a*dy) ;
	}
	bound_var(vz) ;
	bound_var(Bz) ;
}
