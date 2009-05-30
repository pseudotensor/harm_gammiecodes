
#include "decs.h"

void step_ie()
{
	REAL dv ;
	int i,j,k ;
	void bound_var(REAL (*var)[NY+5][NZ+5], REAL tcurr, int sym) ;

	/* update internal energy */
	LOOP {
		dv = (vx[i+1][j][k] - vx[i][j][k])/dx +
		     (vy[i][j+1][k] - vy[i][j][k])/dy +
		     (vz[i][j][k+1] - vz[i][j][k])/dz;

		e[i][j][k] *= (1. - 0.5*dt*(gam-1.)*dv)/(1. + 0.5*dt*(gam-1.)*dv) ;
	}
	bound_var(e,t,SYMMETRIC) ;
}
