
#include "decs.h"

void step_res()
{
	static double *bres ;
	double jx,jy,jz ;
	double **bresx,**bresy,**bresz ;
	int i,j ;
	void bound_var() ;

	/* update internal energy */
	if(fabs(gam - 1.) > 1.e-6) {
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) {
			if(fabs(gam-1.) > 1.e-6) cs = sqrt(gam*(gam-1.)*e[i][j]/rho[i][j]) ;

			jx = (Bz[i][j+1] - Bz[i][j-1])/(2.*a*dy) ;
			jy = (Bz[i+1][j] - Bz[i-1][j])/(2.*a*dx) ;
			jz = (By[i+1][j]+By[i+1][j+1]-By[i-1][j]-By[i-1][j+1])/(4.*a*dx)
				+ (Bx[i][j+1]+Bx[i+1][j+1]-Bx[i][j-1]-Bx[i+1][j-1])/(4.*a*dy) ;

			e[i][j] += dt*(res*cs*dx)*(jx*jx + jy*jy + jz*jz) ;
		}
	}
	bound_var(e) ;

	/* then do magnetic fields */
	bresx = work1 ;
	bresy = work2 ;
	bresz = work3 ;

	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		if(fabs(gam-1.) > 1.e-6) cs = sqrt(gam*(gam-1.)*e[i][j]/rho[i][j]) ;

                bresx[i][j] = dt*(res*cs*dx)*(
			(Bx[i+1][j] - 2.*Bx[i][j] + Bx[i-1][j])/(a*dx*a*dx) +
			(Bx[i][j+1] - 2.*Bx[i][j] + Bx[i][j-1])/(a*dy*a*dy)) ;
                bresy[i][j] = dt*(res*cs*dx)*(
			(By[i+1][j] - 2.*By[i][j] + By[i-1][j])/(a*dx*a*dx) +
			(By[i][j+1] - 2.*By[i][j] + By[i][j-1])/(a*dy*a*dy)) ;
                bresz[i][j] = dt*(res*cs*dx)*(
			(Bz[i+1][j] - 2.*Bz[i][j] + Bz[i-1][j])/(a*dx*a*dx) +
			(Bz[i][j+1] - 2.*Bz[i][j] + Bz[i][j-1])/(a*dy*a*dy)) ;
	}
	bound_var(bresx) ;
	bound_var(bresy) ;
	bound_var(bresz) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		Bx[i][j] += bresx[i][j] ;
		By[i][j] += bresy[i][j] ;
		Bz[i][j] += bresz[i][j] ;
	}
	bound_var(Bx) ;
	bound_var(By) ;
	bound_var(Bz) ;
}
