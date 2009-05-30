
#include "decs.h"

/* driving power spectrum cut off at kmax */
#define KMAX	16

void step_drive()
{
	REAL tmp,rhoa,X,Z ;
	int i,j,k ;

	LOOP {
		X = (i + 0.5)*dx - 0.5*Lx ;
		Z = (k + 0.5)*dz - 0.5*Lz ;

              	tmp = drive*cos(W*t)*Z ;
                rhoa = 0.5*(rho[i][j][k] + rho[i-1][j][k]) ;
                work += dx*dy*dz*rhoa*tmp*vx[i][j][k] ;
                workp += dx*dy*dz*rhoa*drive*sin(W*t)*Z*vx[i][j][k] ;

                vx[i][j][k] += dt*tmp ;
	}

}
