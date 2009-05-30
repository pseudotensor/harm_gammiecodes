
/* enforce periodic boundary conditions */

#include "decs.h"

void bound_var(double **var) 
{
	int i,j ;

	/* y-boundaries */
	for(i=0;i<NX;i++) {
		var[i][-1] = var[i][NY-1] ;
		var[i][NY] = var[i][0] ;
	}
	/* x-boundaries */
	for(j=0;j<NY;j++) {
		var[-1][j] = var[NX-1][j] ;
		var[NX][j] = var[0][j] ;
	}

	/* corners */
	var[-1][-1] = var[NX-1][NY-1] ;
	var[-1][NY] = var[NX-1][0] ;
	var[NX][-1] = var[0][NY-1] ;
	var[NX][NY] = var[0][0] ;
}

