
#include "decs.h"

void sweepym()
{
	REAL X,XP,del,sh,f ;
	int i,j,k,idel,fcase ;
	void bxremap(int i, int k, int idel, int fcase, REAL (*wx)[NX]) ;
	void byremap(int i, int k, int idel, int fcase, REAL (*wy)[NX]) ;
	void bzremap(int i, int k, int idel, int fcase, REAL (*wz)[NX]) ;
	void weights(int i, int fcase, REAL f, REAL sh, REAL (*wx)[NX], REAL (*wy)[NX], REAL (*wz)[NX]) ;
	void bound_var(REAL (*var)[NY+5][NZ+5], REAL tcurr, int sym) ;
	REAL (*bx_dqy)[NY+4][NZ+4],(*bz_dqx)[NY+4][NZ+4],(*bz_dqy)[NY+4][NZ+4] ;
	static REAL wxas[2][NX],wyas[17][NX],wzas[9][NX] ;
	static REAL (*wx)[NX],(*wy)[NX],(*wz)[NX] ;

	wx = (REAL(*)[NX])(&(wxas[0][0])) ;
	wy = (REAL(*)[NX])(&(wyas[0][0])) ;
	wz = (REAL(*)[NX])(&(wzas[0][0])) ;
#if BVANLEER
	bx_dqy = work1 ;
	bz_dqx = work2 ;
	bz_dqy = work3 ;
	
	LOOPP(0,1,0,0,0,0) {
		bx_dqy[i][j][k] = slope_lim(bx[i][j-1][k],bx[i][j][k],bx[i][j+1][k]) ;
	}
	LOOPP(0,0,0,0,0,1) {
		bz_dqx[i][j][k] = slope_lim(bz[i-1][j][k],bz[i][j][k],bz[i+1][j][k]) ;
		bz_dqy[i][j][k] = slope_lim(bz[i][j-1][k],bz[i][j][k],bz[i][j+1][k]) ;
	}
#endif
	sh = q*W*(dx/dy)*dt ;

	/* shift and evolve b_y, preserving div.b = 0  */
	/* must do this first, because byremap assumes that
	   bx,bz still take on their old values */
	for(i=0;i<NX;i++) {
		X = (i + 0.5)*dx - 0.5*Lx ;
		del = -q*W*X*dt/dy ; /* cell shift = (length shift)*(cells/length) or Delta cells = (Delta y)/dy = v0*(Delta t)/dy */
		idel = my_nint(del) ;
		f = del - idel ;
		if(fabs(f) <= sh/2.) fcase = 1 ;
		else if(f > sh/2.) fcase = 2 ;
		else fcase = 3 ;
		weights(i,fcase,f,sh,wx,wy,wz) ; /* calculate weights for all three remaps here */
		for(k=0;k<NZ;k++) byremap(i,k,idel,fcase,wy) ;
	}
	
#if FLUX_REMAP
	for(i=-1;i<NX;i++) {
		XP = i*dx - 0.5*Lx ;
		del = -q*W*XP*dt/dy ;
		/* an updated density is needed to convert back to primitive variables */
		for(k=0;k<NZ;k++) remap_flux(rhotmp,i,k,0,del) ;
	}
#endif
	for(i=0;i<NX;i++) {
		XP = i*dx - 0.5*Lx ;
		del = -q*W*XP*dt/dy ;
	for(k=0;k<NZ;k++) {
#if LINEAR_REMAP
		remap_linear(vx,i,k,-del) ;
#elif FOURIER_REMAP
		remap_fourier(vx,i,k,del) ;
#elif FLUX_REMAP
		remap_flux(vx,i,k,0,del) ;
#endif
	}
	}
#if FLUX_REMAP
	for(i=0;i<NX;i++) {
		X = (i + .5)*dx - 0.5*Lx ;
		del = -q*W*X*dt/dy ;
		for(k=0;k<NZ;k++) remap_flux(rhotmp,i,k,0,del) ;
	}
#endif
	for(i=0;i<NX;i++) {
		X = (i + 0.5)*dx - 0.5*Lx ;
		del = -q*W*X*dt/dy ;
		idel = my_nint(del) ;
		f = del - idel ;
		if(fabs(f) <= sh/2.) fcase = 1 ;
		else if(f > sh/2.) fcase = 2 ;
		else fcase = 3 ;
	for(k=0;k<NZ;k++) {
		bxremap(i,k,idel,fcase,wx) ;
		bzremap(i,k,idel,fcase,wz) ;
#if LINEAR_REMAP
		remap_linear(rho,i,k,-del) ;
		remap_linear(e,i,k,-del) ;
		remap_linear(vy,i,k,-del) ;
		remap_linear(vz,i,k,-del) ;
#elif FOURIER_REMAP
		remap_fourier(rho,i,k,del) ;
		remap_fourier(e,i,k,del) ;
		remap_fourier(vy,i,k,del) ;
		remap_fourier(vz,i,k,del) ;
#elif FLUX_REMAP
		remap_flux(e,i,k,0,del) ;
		remap_flux(vy,i,k,0,del) ;
		remap_flux(vz,i,k,0,del) ;
		remap_flux(rho,i,k,0,del) ; /* need to update rho last */
#endif
	}
	}

	bound_var(rho,t+dt,SYMMETRIC) ; /* need to bound rho first */
	bound_var(e,t+dt,SYMMETRIC) ;
	bound_var(vx,t+dt,SYMMETRIC) ;
	bound_var(vy,t+dt,SYMMETRIC) ;
	bound_var(vz,t+dt,ANTISYMMETRIC) ;
	bound_var(by,t+dt,SYMMETRIC) ;
	bound_var(bz,t+dt,ANTISYMMETRIC) ;
        bound_var(bx,t+dt,SYMMETRIC) ;
}
