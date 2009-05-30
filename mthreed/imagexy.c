/*
	produces an "r8" file.  
*/

#include "decs.h"

void image(FILE *fp)
{
	float iq,liq,a,b,lmax,lmin,X ;
	int i,j ;
	static int slice = NZ/2 ;

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
#elif BY0
	lmax = 0. ;
	lmin = -q*W*B0*tf ;
#endif
	a = 256./(lmax - lmin) ;
	b = -a*lmin ;

	//for(j=NY+1;j>=-2;j--) 
	//for(i=-2;i<=NX+1;i++) {
	for(j=NY-1;j>=0;j--)
	for(i=0;i<=NX-1;i++) {
		X = (i + .5)*dx - .5 ;
#if LINEAR
		iq = vz[i][j][slice] ;
		iq = vx[i][j][slice] ;
		iq = vy[i][j][slice] ;
		//iq = rho[i][j][slice] - 1. ;
#elif RHOCOL || STRIPE
		iq = rho[i][j][slice] ;
#elif FLOOP
		iq = (by[i][j][slice] - by[i-1][j][slice])/dx - (bx[i][j][slice] - bx[i][j-1][slice])/dy ;
#elif BY0
		iq = by[i][j][slice] ;
#endif
		liq = a*iq + b ;
		if(liq > 255.) liq = 255. ;
		if(liq < 0.) liq = 0. ;
		fprintf(fp,"%c",(char)((int)liq)) ;
	}
}
