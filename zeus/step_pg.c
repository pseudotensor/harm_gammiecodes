
#include "decs.h"

/* pressure and gravity step */
void step_pg()
{
	double rhoa,bxa,bxp,bxm,bya,byp,bym,bza,bzp,bzm ;
	int i,j ;
	void bound_var() ;

	/* first pressure gradient */
	if(fabs(gam - 1.) > 1.e-6) {
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) {
			rhoa = 0.5*(rho[i][j]+rho[i-1][j]) ;
			vx[i][j] += -(dt/(a*dx*rhoa))*(gam-1.)*(e[i][j]-e[i-1][j]) ;
			rhoa = 0.5*(rho[i][j]+rho[i][j-1]) ;
			vy[i][j] += -(dt/(a*dy*rhoa))*(gam-1.)*(e[i][j]-e[i][j-1]) ;
		}
	}
	else {
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) {
			rhoa = 0.5*(rho[i][j]+rho[i-1][j]) ;
			vx[i][j] += -(dt/(a*dx*rhoa))*cs*cs*(rho[i][j]-rho[i-1][j]) ;
			rhoa = 0.5*(rho[i][j]+rho[i][j-1]) ;
			vy[i][j] += -(dt/(a*dy*rhoa))*cs*cs*(rho[i][j]-rho[i][j-1]) ;
		}
	}

	/* gravitational acceleration */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vx[i][j] += -(dt/(a*dx))*(pot[i][j]-pot[i-1][j]) ;
		vy[i][j] += -(dt/(a*dy))*(pot[i][j]-pot[i][j-1]) ;
	}

	/* magnetic pressure gradient */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		rhoa = 0.5*(rho[i][j]+rho[i-1][j]) ;

		bya = 0.25*(By[i][j]+By[i-1][j]+By[i-1][j+1]+By[i][j+1]) ;
		byp = 0.5*(By[i][j]+By[i][j+1]) ;
		bym = 0.5*(By[i-1][j]+By[i-1][j+1]) ;

		bza = 0.5*(Bz[i][j]+Bz[i-1][j]) ;
		bzp = Bz[i][j] ;
		bzm = Bz[i-1][j] ;

		vx[i][j] += -(dt/(a*dx*rhoa))*(
			bya*(byp - bym) +
			bza*(bzp - bzm)
			) ;
			
		rhoa = 0.5*(rho[i][j]+rho[i][j-1]) ;

		bxa = 0.25*(Bx[i][j]+Bx[i][j-1]+Bx[i+1][j-1]+Bx[i+1][j]) ;
		bxp = 0.5*(Bx[i][j]+Bx[i+1][j]) ;
		bxm = 0.5*(Bx[i][j-1]+Bx[i+1][j-1]) ;

		bza = 0.5*(Bz[i][j]+Bz[i][j-1]) ;
		bzp = Bz[i][j] ;
		bzm = Bz[i][j-1] ;

		vy[i][j] += -(dt/(a*dy*rhoa))*(
			bxa*(bxp - bxm) +
			bza*(bzp - bzm)
			) ;
	}

	/* only needs to be done once */
	bound_var(vx) ;
	bound_var(vy) ;

}
