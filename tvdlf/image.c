/* 
	produces an "r8" file.  
*/

#include "decs.h"

void image(FILE *fp)
{
	double iq,liq,a,b,lmax,lmin ;
	int i,j ;

	/* density mapping is logarithmic, in 255 steps
	   between e^lmax and e^lmin */
	lmax = 2. ;
	lmin = 0. ;
	a = 256./(lmax - lmin) ;
	b = -a*lmin ;

	for(j=NY-1;j>=0;j--) 
	for(i=0;i<NX;i++) {
		iq = (fabs(p[i][j][UU]) + SMALL)*(gam - 1.)  ;
		liq = a*log(iq) + b ;

		if(liq > 255.) liq = 255. ;
		if(liq < 0.) liq = 0. ;
		fprintf(fp,"%c",(char)((int)liq)) ;
	}
}

