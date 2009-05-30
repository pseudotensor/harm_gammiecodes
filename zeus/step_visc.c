
#include "decs.h"

void step_visc()
{
	static double **visx,**visy ;
	double dvx,dvy,qlx,qvnrx,qly,qvnry ;
	double rhoa,rhob ;
	int i,j ;
	void bound_var() ;

	visx = work8 ;
	visy = work9 ;
	
	/* find viscous stresses */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		/* del v, at zone center */
		dvx = vx[i+1][j] - vx[i][j] ;
		dvy = vy[i][j+1] - vy[i][j] ;
		
		/* linear viscosity */
		if(fabs(gam - 1.) > 1.e-6) cs = sqrt(gam*(gam-1.)*e[i][j]/rho[i][j]) ;
		qlx = -nu_l*rho[i][j]*cs*dvx ;
		qly = -nu_l*rho[i][j]*cs*dvy ;

		/* von neumann,richtmyer viscosity */
		if(dvx  < 0) {
			qvnrx = nu_vnr*rho[i][j]*dvx*dvx ;
		}
		else qvnrx = 0. ;
		if(dvy  < 0) {
			qvnry = nu_vnr*rho[i][j]*dvy*dvy ;
		}
		else qvnry = 0. ;

		visx[i][j] = qlx + qvnrx ;
		visy[i][j] = qly + qvnry ;
	}
	bound_var(visx) ;
	bound_var(visy) ;

	/* update velocity, internal energy */
	if(fabs(gam - 1.) > 1.e-6) {
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) {
			e[i][j] += 
				-dt*visx[i][j]*(vx[i+1][j] - vx[i][j])/(a*dx) 
				- dt*visy[i][j]*(vy[i][j+1] - vy[i][j])/(a*dy) ;
		}
	}
	bound_var(e) ;

	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		rhoa = 0.5*(rho[i][j]+rho[i-1][j]) ;
		vx[i][j] += -dt*(visx[i][j]-visx[i-1][j])/(a*dx*rhoa) ;
		rhob = 0.5*(rho[i][j]+rho[i][j-1]) ;
		vy[i][j] += -dt*(visy[i][j]-visy[i][j-1])/(a*dy*rhob) ;
	}
	bound_var(vx) ;
	bound_var(vy) ;
}
