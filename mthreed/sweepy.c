
#include "decs.h"

void sweepy()
{
	REAL vpp,mdotp ;
	int i,j,k ;
	REAL (*px)[NY+4][NZ+4], (*py)[NY+4][NZ+4], (*pz)[NY+4][NZ+4] ;
	REAL (*dq)[NY+4][NZ+4], (*mdot)[NY+4][NZ+4], (*fl)[NY+4][NZ+4] ;
	void bound_var(REAL (*var)[NY+5][NZ+5], REAL tcurr, int sym) ;

	px = work1 ;
	py = work2 ;
	pz = work3 ;
	dq = work4 ;
	mdot = work5 ;
	fl = work6 ;

	/* transform to momenta */
	LOOP {
		px[i][j][k] = 0.5*(rho[i-1][j][k]+rho[i][j][k])*vx[i][j][k] ;
		py[i][j][k] = 0.5*(rho[i][j-1][k]+rho[i][j][k])*vy[i][j][k] ;
		pz[i][j][k] = 0.5*(rho[i][j][k-1]+rho[i][j][k])*vz[i][j][k] ;
	}

	/* create normalized transport variable */
	LOOPP(2,2,2,2,2,2) vy[i][j][k] *= dt/dy ;

	/* first do mass */
	LOOPP(1,0,2,1,1,0) {
		dq[i][j][k] = slope_lim(rho[i][j-1][k],rho[i][j][k],rho[i][j+1][k]) ;
	}
	LOOPP(1,0,1,1,1,0) {
		mdot[i][j][k] = (vy[i][j][k] > 0.) ?
			vy[i][j][k]*(rho[i][j-1][k] + (1. - vy[i][j][k])*dq[i][j-1][k]) :
			vy[i][j][k]*(rho[i][j][k] - (1. + vy[i][j][k])*dq[i][j][k]) ;

	}

	/* then energy */
	LOOPP(0,0,1,1,0,0) {
		dq[i][j][k] = slope_lim(e[i][j-1][k]/rho[i][j-1][k],
					e[i][j][k]/rho[i][j][k],
					e[i][j+1][k]/rho[i][j+1][k]) ;
	}
	LOOPP(0,0,0,1,0,0) {
		fl[i][j][k] = (vy[i][j][k] > 0.) ?
			mdot[i][j][k]*(e[i][j-1][k]/rho[i][j-1][k] + (1. - vy[i][j][k])*dq[i][j-1][k]) :
			mdot[i][j][k]*(e[i][j][k]/rho[i][j][k] - (1. + vy[i][j][k])*dq[i][j][k]) ;
	}
	LOOP e[i][j][k] += -(fl[i][j+1][k]-fl[i][j][k]) ;

	/* vx: because of different centering, indices differ */
	LOOPP(0,0,1,1,0,0) {
		dq[i][j][k] = slope_lim(vx[i][j-1][k],vx[i][j][k],vx[i][j+1][k]) ;
	}
	LOOPP(0,0,0,1,0,0) {
		vpp = 0.5*(vy[i][j][k]+vy[i-1][j][k]) ;
		mdotp = 0.5*(mdot[i][j][k]+mdot[i-1][j][k]) ;
		fl[i][j][k] = (vpp > 0.) ?
			mdotp*(vx[i][j-1][k] + (1. - vpp)*dq[i][j-1][k]) :
			mdotp*(vx[i][j][k] - (1. + vpp)*dq[i][j][k]) ;
	}
	LOOP px[i][j][k] += -(fl[i][j+1][k]-fl[i][j][k]) ;

	/* vy */
	LOOPP(0,0,1,1,0,0) {
		dq[i][j][k] = slope_lim(vy[i][j-1][k],vy[i][j][k],vy[i][j+1][k]) ;
	}
	LOOPP(0,0,1,0,0,0) {
		vpp = 0.5*(vy[i][j][k]+vy[i][j+1][k]) ;
		mdotp = 0.5*(mdot[i][j][k]+mdot[i][j+1][k]) ;
		fl[i][j][k] = (vpp > 0.) ?
			mdotp*(vy[i][j][k] + (1. - vpp)*dq[i][j][k]) :
			mdotp*(vy[i][j+1][k] - (1. + vpp)*dq[i][j+1][k]) ;
	}
	LOOP py[i][j][k] += -(fl[i][j][k]-fl[i][j-1][k])*(dy/dt) ;

	/* vz */
	LOOPP(0,0,1,1,0,0) {
		dq[i][j][k] = slope_lim(vz[i][j-1][k],vz[i][j][k],vz[i][j+1][k]) ;
	}
	LOOPP(0,0,0,1,0,0) {
		vpp = 0.5*(vy[i][j][k]+vy[i][j][k-1]) ;
		mdotp = 0.5*(mdot[i][j][k]+mdot[i][j][k-1]) ;
		fl[i][j][k] = (vpp > 0.) ?
			mdotp*(vz[i][j-1][k] + (1. - vpp)*dq[i][j-1][k]) :
			mdotp*(vz[i][j][k] - (1. + vpp)*dq[i][j][k]) ;
	}
	LOOP pz[i][j][k] += -(fl[i][j+1][k]-fl[i][j][k]) ;

	LOOP {
		rho[i][j][k] += -(mdot[i][j+1][k]-mdot[i][j][k]) ;
		if(rho[i][j][k] < MIN) {
			fprintf(stderr,"y small density: %d %g\n",i,rho[i][j][k]) ;
		}
	}
	bound_var(rho,t,SYMMETRIC) ;

	/* return to velocities */
	LOOP {
		vx[i][j][k] = 2.*px[i][j][k]/(rho[i][j][k]+rho[i-1][j][k]) ;
		vy[i][j][k] = 2.*py[i][j][k]/(rho[i][j][k]+rho[i][j-1][k]) ;
		vz[i][j][k] = 2.*pz[i][j][k]/(rho[i][j][k]+rho[i][j][k-1]) ;
	}
	bound_var(vx,t,SYMMETRIC) ;
	bound_var(vy,t,SYMMETRIC) ;
	bound_var(vz,t,ANTISYMMETRIC) ;
	bound_var(e,t,SYMMETRIC) ;
}
