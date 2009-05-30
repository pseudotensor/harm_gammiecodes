/*
	produces an "r8" file.  
*/

#include "decs.h"

void image(FILE *fp)
{
	float iq,liq,a,b,lmax,lmin,X ;
	int i,j,k ;
	static int slice = NX/2 ;

	/* density mapping is logarithmic, in 255 steps
	   between e^lmax and e^lmin */
#if LINEAR
	lmax = pert ;
	lmin = -pert ;
#elif RHOCOL || STRIPE
	lmax = 3. ;
	lmin = 1. ;
#elif FLOOP
	lmax = 0.08 ;
	lmin = -0.04 ;
#endif
	a = 256./(lmax - lmin) ;
	b = -a*lmin ;

	for(k=NZ-1;k>=0;k--)
	for(j=0;j<=NY-1;j++) {
		X = (i + .5)*dx - .5 ;
#if LINEAR
		iq = vz[slice][j][k] ;
		iq = vx[slice][j][k] ;
		iq = vy[slice][j][k] ;
#elif RHOCOL || STRIPE
		iq = rho[slice][j][k] ;
#elif FLOOP
		iq = -(bx[slice][j][k] - bx[slice][j-1][k])/dy ;
#endif
		liq = a*iq + b ;
		if(liq > 255.) liq = 255. ;
		if(liq < 0.) liq = 0. ;
		fprintf(fp,"%c",(char)((int)liq)) ;
	}
}
