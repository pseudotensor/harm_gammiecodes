
#include "decs.h"

/* special diagnostics subroutine */

void spec_diag()
{
	static double px0,py0,pz0 ;
	double px,py,pz ;
	double mono ;
	int i,j ;
	static int firstc = 1;

	/* mean momenta */
	px = py = pz = 0. ;
	mono = 0. ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		px += rho[i][j]*0.5*(vx[i][j] + vx[i+1][j]) ;
		py += rho[i][j]*0.5*(vy[i][j] + vy[i][j+1]) ;
		pz += rho[i][j]*vz[i][j] ;
		mono += (Bx[i+1][j]-Bx[i][j])/dx +
			(By[i][j+1]-By[i][j])/dy ;
	}
	px /= NX*NY ;
	py /= NX*NY ;
	pz /= NX*NY ;
	mono /= NX*NY ;

	fprintf(stderr,"%10.5g %10.5g %10.5g %10.5g\n",
		px-px0,py-py0,pz-pz0,mono) ;

	if(firstc) {
		firstc = 0 ;
		px0 = px ;
		py0 = py ;
		pz0 = pz ;
	}

}
