

/* magnetic pressure step */

#include "decs.h"

void step_pb()
{
	double rhoa,bza,bya,bxa,byp,bym,bxp,bxm ;
	int i,j ;
	void bound_var() ;

	/* magnetic pressure gradient */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		bya = 0.25*(By[i-1][j] + By[i][j] + By[i][j+1] + By[i-1][j+1]) ;
		byp = 0.5*(By[i][j] + By[i][j+1]) ;
		bym = 0.5*(By[i-1][j] + By[i-1][j+1]) ;
		bza = 0.5*(Bz[i][j]+Bz[i-1][j]) ;
		rhoa = 0.5*(rho[i][j]+rho[i-1][j]) ;
		vx[i][j] += -(dt/(a*dx*rhoa))*(
			bya*(byp - bym) +
			bza*(Bz[i][j] - Bz[i-1][j])
			) ;

		bxa = 0.25*(Bx[i][j-1] + Bx[i+1][j-1] + Bx[i+1][j] + Bx[i][j]) ;
		bxp = 0.5*(Bx[i][j] + Bx[i+1][j]) ;
		bxm = 0.5*(Bx[i][j-1] + Bx[i+1][j-1]) ;
		bza = 0.5*(Bz[i][j]+Bz[i][j-1]) ;
		rhoa = 0.5*(rho[i][j]+rho[i][j-1]) ;
		vy[i][j] += -(dt/(a*dy*rhoa))*(
			bxa*(bxp - bxm) +
			bza*(By[i][j] - Bz[i][j-1])
			) ;
	}
	bound_var(vx) ;
	bound_var(vy) ;

}
