
#include "decs.h"

void step_ie()
{
	double dv ;
	int i,j ;
	void bound_var() ;

	/* update internal energy */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		dv = (vx[i+1][j] - vx[i][j])/(a*dx) +
		     (vy[i][j+1] - vy[i][j])/(a*dy) ;

		e[i][j] *= (1. - 0.5*dt*(gam-1.)*dv)/(1. + 0.5*dt*(gam-1.)*dv) ;
	}
	bound_var(e) ;
}
