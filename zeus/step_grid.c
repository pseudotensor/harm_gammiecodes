
/* add source terms due to grid expansion and contraction */

#include "decs.h"

void step_grid()
{
	int i,j ;
	double dv ;
	void bound_var() ;

	/* magnetic fields */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		Bx[i][j] += -2.*dt*da*Bx[i][j]/a ;
		By[i][j] += -2.*dt*da*By[i][j]/a ;
		Bz[i][j] += -2.*dt*da*Bz[i][j]/a ;
	}
	bound_var(Bx) ;
	bound_var(By) ;
	bound_var(Bz) ;

	/* density */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		rho[i][j] += -3.*dt*da*rho[i][j]/a ;
	}
	bound_var(rho) ;

	/* internal energy */
	if(fabs(gam - 1.) > 1.e-6) {
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) {
			dv = 3.*da/a ;
			e[i][j] *= (1 - 0.5*dt*(gam-1.)*dv)/(1 + 0.5*dt*(gam-1.)*dv) ;
		}
	}
	bound_var(e) ;

	/* velocities */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vx[i][j] += -dt*da*vx[i][j]/a ;
		vy[i][j] += -dt*da*vy[i][j]/a ;
		vz[i][j] += -dt*da*vz[i][j]/a ;
	}
	bound_var(vx) ;
	bound_var(vy) ;
	bound_var(vz) ;
}
