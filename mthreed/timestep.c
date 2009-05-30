
#include "decs.h"

void timestep() 
{
	int i,j,k,imax,jmax,kmax ;
	REAL u,va,X ;
	REAL idt1,idt2,idt3,idt4,idt5,idt6,idt7,idt8,idt9 ;
	REAL dtinv ;
	REAL dtinv_max,dt3m,dt4m ;
	REAL dv,d,vxpot,vypot,dvy,fast ;
	static REAL dtlast ;
	void dump() ;

	if(dx <= dy && dx <= dz) d = dx ;
	else if (dy <= dx && dy <= dz) d = dy ;
	else d = dz ;

        /* SSBC timestep criterion */
	idt2 = 2.*q*W ;
	idt2 *= idt2 ;
	idt2 = 0. ;

	dtinv_max = 0. ;
	dtlast = dt ;
	LOOP {
		if(rho[i][j][k] < 0) {
			fprintf(stderr,"negative density error\n") ;
			dump(stdout) ;
			exit(4) ;
		}
		
		if(fabs(gam - 1.) > 1.e-6) cs = sqrt(gam*(gam-1.)*e[i][j][k]/rho[i][j][k]) ;

		/* fast speed; not centered */
		fast = cs*cs +
			(bx[i][j][k]*bx[i][j][k] + by[i][j][k]*by[i][j][k] + bz[i][j][k]*bz[i][j][k])/rho[i][j][k]  ;
		idt1 = fast/(d*d) ;

		/* x-velocity
		X = (i + 0.5)*dx - .5*Lx ;
		vxpot = vx[i][j][k]*vx[i][j][k] + 4.*vy[i][j][k]*vy[i][j][k] ;
		idt3 = vxpot/(dx*dx) ;

                /* y-velocity; term to eliminate epicyclic bias
		vypot = vxpot/4. ;
		idt4 = vypot/(dy*dy) ; */

		/* x-velocity */
		idt3 = vx[i][j][k]/dx ;
		idt3 *= idt3 ;

		/* y-velocity */
		idt4 = vy[i][j][k]/dy ;
		idt4 *= idt4 ;

		/* z-velocity */
		idt5 = vz[i][j][k]/dz ;
		idt5 *= idt5 ;

		/* linear viscosity */
		idt6 = 4.*nu_l*cs/d ;
		idt6 *= idt6 ;

		/* VNR viscosity */
		dv = vx[i+1][j][k]-vx[i][j][k] ;
		idt7 = 4.*nu_vnr*dv/dx ;
		idt7 *= idt7 ;

		/* VNR viscosity */
		dv = vy[i][j+1][k]-vy[i][j][k] ;
		idt8 = 4.*nu_vnr*dv/dy ;
		idt8 *= idt8 ;

		/* VNR viscosity */
		//dv = vy[i][j][k+1]-vy[i][j][k] ;
		dv = vz[i][j][k+1]-vz[i][j][k] ;
		idt9 = 4.*nu_vnr*dv/dz ;
		idt9 *= idt9 ;

		dtinv = idt1 + idt2 + idt3 + idt4 + idt5 + idt6 + idt7 + idt8 + idt9 ;

		if(dtinv > dtinv_max) {
			dtinv_max = dtinv ;
			imax = i ;
			jmax = j ;
			kmax = k ;
			dt3m = idt3 ;
			dt4m = idt4 ;
		}
	}
	/*
	fprintf(stderr,"%d %d %d %g %g\n",imax,jmax,kmax,dt3m,dt4m) ;
	fprintf(stderr,"\n") ;
	*/

	dt = cour/sqrt(dtinv_max) ;

	/* don't increase timestep by too much */
	if(dt > 1.3*dtlast) dt = 1.3*dtlast ;

	/* don't step beyond end of run */
	if(t + dt > tf) dt = tf - t ;
}

