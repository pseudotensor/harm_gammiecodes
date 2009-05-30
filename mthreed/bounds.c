
/* enforce boundary conditions */

#define	REFLECTING	0
#define PERIODIC	1

#include "decs.h"

void bound_var(REAL (*var)[NY+5][NZ+5], REAL tcurr, int sym)
{
	int i,j,k,lo ;
	REAL del ;
	
	if(var==rho) lo = -3 ;
	else lo = -2 ;

	/* x-boundaries: shearing box BCs  */
#if FLUX_REMAP
#else
	for(j=0;j<NY;j++)
	for(k=0;k<NZ;k++) {
		for(i=lo;i<=-1;i++) var[i][j][k] = var[i+NX][j][k] ;
		var[NX][j][k] = var[0][j][k] ;
		var[NX+1][j][k] = var[1][j][k] ;
	}
#endif
	if(q!=0.) {
	del = -q*W*Lx*tcurr/dy ;
	for(k=0;k<NZ;k++) {
#if LINEAR_REMAP
		for(i=lo;i<=-1;i++) remap_linear(var,i,k,del) ;
		remap_linear(var,NX,k,-del) ;
		remap_linear(var,NX+1,k,-del) ;
#elif FOURIER_REMAP
		for(i=lo;i<=-1;i++) remap_fourier(var,i,k,del) ;
		remap_fourier(var,NX,k,del) ;
		remap_fourier(var,NX+1,k,del) ;
#elif FLUX_REMAP
		if(var==vx||var==vy||var==bx||var==by) remap_flux(rhotmp,-1,k,NX,-del) ;
		if(var==vx||var==bx) remap_flux(rhotmp,-2,k,NX,-del) ;
		remap_flux(var,-1,k,NX,-del) ;
		if(var==vx||var==vy||var==bx||var==by) remap_flux(rhotmp,-2,k,NX,-del) ;
		if(var==vx||var==bx) remap_flux(rhotmp,-3,k,NX,-del) ;
		remap_flux(var,-2,k,NX,-del) ;
		if(var==rho) remap_flux(var,-3,k,NX,-del) ;
		if(var==vx||var==vy||var==bx||var==by) remap_flux(rhotmp,NX,k,-NX,del) ;
        	if(var==vx) remap_flux(rhotmp,NX-1,k,-NX,del) ;
		remap_flux(var,NX,k,-NX,del) ;
		if(var==vx||var==vy||var==bx||var==by) remap_flux(rhotmp,NX+1,k,-NX,del) ;
        	if(var==vx||var==bx) remap_flux(rhotmp,NX,k,-NX,del) ;
		remap_flux(var,NX+1,k,-NX,del) ;
#endif
	}
	}
	
	if(var==bx) {
	for(j=0;j<NY;j++)
	for(k=0;k<NZ;k++) {
		bx[NX][j][k] = bx[NX-1][j][k] + dx*(
			    + (by[NX-1][j][k] - by[NX-1][j+1][k])/dy
			    + (bz[NX-1][j][k] - bz[NX-1][j][k+1])/dz) ;
	}
	}
	
#if SWEEPYM
#else
        if(var==vy&&q!=0.) for(j=0;j<NY;j++)
	for(k=0;k<NZ;k++) {
                var[NX][j][k] += -q*W*Lx ;
                var[NX+1][j][k] += -q*W*Lx ;
                var[-1][j][k] += q*W*Lx ;
                var[-2][j][k] += q*W*Lx ;
        }
#endif

	/* y-boundaries: periodic */
	for(i=lo;i<NX+2;i++)
	for(k=0;k<NZ;k++) {
		for(j=lo;j<=-1;j++) var[i][j][k] = var[i][j+NY][k] ;
		var[i][NY][k] = var[i][0][k] ;
		var[i][NY+1][k] = var[i][1][k] ;
	}

#if PERIODIC
	for(i=lo;i<NX+2;i++)
	for(j=lo;j<NY+2;j++) {
		for(k=lo;k<=-1;k++) var[i][j][k] = var[i][j][k+NZ] ;
		var[i][j][NZ] = var[i][j][0] ;
		var[i][j][NZ+1] = var[i][j][1] ;
	}
#endif

#if REFLECTING
	/* z-boundaries: reflecting, for symmetric variables rho,vx,vy */
	if(sym == SYMMETRIC) {
		for(i=lo;i<NX+2;i++)
		for(j=lo;j<NY+2;j++) {
			var[i][j][-1] = var[i][j][0] ;
			var[i][j][-2] = var[i][j][1] ;
			if(var==rho) var[i][j][-3] = var[i][j][2] ;
			var[i][j][NZ] = var[i][j][NZ-1] ;
			var[i][j][NZ+1] = var[i][j][NZ-2] ;
		}
	}
	else {	/* for antisymmetric variables, e.g. vz */
		for(i=0;i<=NX+1;i++)
		for(j=0;j<=NY+1;j++) {
			var[i][j][-2] = -var[i][j][2] ;
			var[i][j][-1] = -var[i][j][1] ;
			var[i][j][0] = 0. ;
			var[i][j][NZ] = 0. ;
			var[i][j][NZ+1] = -var[i][j][NZ-1] ;
		}
	}
#endif

#if OUTFLOW
#endif
}

