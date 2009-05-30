
#include "decs.h"

void sweepx()
{
	double vpp,mdotp ;
	double **dq,**vp ;
	double **mdot ;
	double **px,**py,**pz ;
	double **u,**Bzr ;
	int i,j ;
	void bound_var() ;
	void dqx_calc(double **var, double **dq) ;

	/* transform to momenta */
	px = work1 ;
	py = work2 ;
	pz = work3 ;
	vp = work9 ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		px[i][j] = 0.5*(rho[i-1][j]+rho[i][j])*vx[i][j] ;
		py[i][j] = 0.5*(rho[i][j-1]+rho[i][j])*vy[i][j] ;
		pz[i][j] = rho[i][j]*vz[i][j] ;

		/* create dimensionless transport variable */
		vp[i][j] = vx[i][j]*dt/(a*dx) ;
	}
	bound_var(vp) ;

	mdot = work4 ;
	dq = work6 ;

	/* first do mass */
	dqx_calc(rho,dq) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		if(vp[i][j] > 0.) 
			fl[i][j] = vp[i][j]*(rho[i-1][j] + (1. - vp[i][j])*dq[i-1][j]) ;
		else	
			fl[i][j] = vp[i][j]*(rho[i][j] - (1. + vp[i][j])*dq[i][j]) ;
		mdot[i][j] = fl[i][j] ;
	}
	bound_var(fl) ;
	bound_var(mdot) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		rho[i][j] += -(fl[i+1][j] - fl[i][j]) ;
		if(rho[i][j] < MIN) {
			fprintf(stderr,"small density: %d %g\n",i,rho[i][j]) ;
		}
	}
	bound_var(rho) ;

	/* then energy */
	u = work7 ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		u[i][j] = e[i][j]/rho[i][j] ;
	}
	bound_var(u) ;
	dqx_calc(u,dq) ;

	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		if(vp[i][j] > 0.) 
			fl[i][j] = mdot[i][j]*(u[i-1][j] + 
				(1. - vp[i][j])*dq[i-1][j]) ;
		else	
			fl[i][j] = mdot[i][j]*(u[i][j] - 
				(1. + vp[i][j])*dq[i][j]) ;
	}
	bound_var(fl) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		e[i][j] += -(fl[i+1][j]-fl[i][j]) ;
	}
	bound_var(e) ;

	/* magnetic field: bz */
	Bzr = work8 ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		Bzr[i][j] = Bz[i][j]/rho[i][j] ;
	}
	bound_var(Bzr) ;
	dqx_calc(Bzr,dq) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		if(vp[i][j] > 0) 
			fl[i][j] = mdot[i][j]*(Bzr[i-1][j] + 
				(1. - vp[i][j])*dq[i-1][j]) ;
		else	
			fl[i][j] = mdot[i][j]*(Bzr[i][j] - 
				(1. + vp[i][j])*dq[i][j]) ;
	}
	bound_var(fl) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		Bz[i][j] += -(fl[i+1][j]-fl[i][j]) ;
	}
	bound_var(Bz) ;

	/* vx: because of different centering, indices differ */
	dqx_calc(vx,dq) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vpp = 0.5*(vp[i+1][j]+vp[i][j]) ;
		mdotp = 0.5*(mdot[i+1][j]+mdot[i][j]) ;
		if(vpp > 0) 
			fl[i][j] = mdotp*(vx[i][j] + (1. - vpp)*dq[i][j]) ;
		else	
			fl[i][j] = mdotp*(vx[i+1][j] - (1. + vpp)*dq[i+1][j]) ;
	}
	bound_var(fl) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		px[i][j] += -(fl[i][j]-fl[i-1][j]) ;
	}
	bound_var(px) ;

	/* vy */
	dqx_calc(vy,dq) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		mdotp = 0.5*(mdot[i][j]+mdot[i][j-1]) ;
		vpp = 0.5*(vp[i][j]+vp[i][j-1]) ;
		if(vpp > 0.) 
			fl[i][j] = mdotp*(vy[i-1][j] + 
				(1. - vpp)*dq[i-1][j]) ;
		else	
			fl[i][j] = mdotp*(vy[i][j] - 
				(1. + vpp)*dq[i][j]) ;
	}
	bound_var(fl) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		py[i][j] += -(fl[i+1][j]-fl[i][j]) ;
	}
	bound_var(py) ;

	/* vz */
	dqx_calc(vz,dq) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		if(vp[i][j] > 0) 
			fl[i][j] = mdot[i][j]*(vz[i-1][j] + 
				(1. - vp[i][j])*dq[i-1][j]) ;
		else	
			fl[i][j] = mdot[i][j]*(vz[i][j] - 
				(1. + vp[i][j])*dq[i][j]) ;
	}
	bound_var(fl) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		pz[i][j] += -(fl[i+1][j]-fl[i][j]) ;
	}
	bound_var(pz) ;

	/* return to velocities */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vx[i][j] = 2.*px[i][j]/(rho[i][j]+rho[i-1][j]) ;
		vy[i][j] = 2.*py[i][j]/(rho[i][j]+rho[i][j-1]) ;
		vz[i][j] = pz[i][j]/rho[i][j] ;
	}
	bound_var(vx) ;
	bound_var(vy) ;
	bound_var(vz) ;
}

void dqx_calc(double **var,double **dq)
{
	int i,j ;
	double Dqp,Dqm,pr ;

        for(i=0;i<NX;i++) 
        for(j=0;j<NY;j++) {
		Dqp = var[i+1][j] - var[i][j] ;
		Dqm = var[i][j] - var[i-1][j] ;

                pr = Dqp*Dqm ;
		dq[i][j] = (pr > 0.) ? pr/(Dqm + Dqp) : 0. ;
        }
        bound_var(dq) ;
}

