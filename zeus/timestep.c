
#include "decs.h"

void timestep() 
{
	int i,j ;
	double u,va ;
	double idt1,idt2,idt3,idt4,idt5,idt6,idt7,idt8 ;
	double **dtinv,dtinv_max ;
	double bxa,bya,dv,d ;
	static double dtlast ;
	void dump() ;

	if(dx > dy) d = dy ;
	else d = dx ;

	dtinv = work1 ;

	dtinv_max = 0. ;
	dtlast = dt ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		u = e[i][j] ;
		if(u < 0) {
			fprintf(stderr,"negative internal energy error\n") ;
			dump(stdout) ;
			exit(3) ;
		}
		if(rho[i][j] < 0) {
			fprintf(stderr,"negative density error\n") ;
			dump(stdout) ;
			exit(4) ;
		}

		/* sound speed */
		if(fabs(gam - 1.) > SMALL) cs = sqrt(gam*(gam-1.)*u/rho[i][j]) ;
		idt1 = cs/(a*d) ;

		/* x-velocity */
		idt2 = (vx[i][j]+vx[i+1][j])/(2.*a*dx) ;

		/* y-velocity */
		idt3 = (vy[i][j]+vy[i][j+1])/(2.*a*dy) ;

		/* alfven velocity */
		bxa = 0.5*(Bx[i][j] + Bx[i+1][j]) ;
		bya = 0.5*(By[i][j] + By[i][j+1]) ;
		va = sqrt((bxa*bxa + bya*bya + Bz[i][j]*Bz[i][j])/rho[i][j]) ;
		idt4 = va/(a*d) ;

		/* linear viscosity */
		idt5 = 4.*nu_l*cs/(a*d) ;

		/* VNR viscosity */
		dv = vx[i+1][j]-vx[i][j] ;
		idt6 = (4.*nu_vnr*dv)/(a*dx) ;

		/* VNR viscosity */
		dv = vy[i][j+1]-vy[i][j] ;
		idt7 = 4.*nu_vnr*dv/(a*dy) ;

		/* resistivity */
		idt8 = 4.*res*cs/(a*d) ;

		dtinv[i][j] = idt1*idt1 +
			idt2*idt2 +
			idt3*idt3 +
			idt4*idt4 +
			idt5*idt5 +
			idt6*idt6 +
			idt7*idt7 +
			idt8*idt8 ;
	}

	
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		if(dtinv[i][j] > dtinv_max) dtinv_max = dtinv[i][j] ;
	}

	dt = cour/sqrt(dtinv_max) ;

	/* don't increase timestep by too much */
	if(dt > 1.3*dtlast) dt = 1.3*dtlast ;

	/* don't step beyond next driving time */
	if(t + dt > tdrive && ampf > 0.) dt = tdrive - t + SMALL ;

	/* don't step beyond end of run */
	if(t + dt > tf) dt = tf - t ;

}
