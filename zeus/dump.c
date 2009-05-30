
#include "decs.h"

#define MOVIE	0
#define NO_MOV	1

void dump(fp)
FILE *fp ;
{
	int i,j ;
	double vxa,vya,vza ;
	double Bxa,Bya,Bza ;
	double eb,ek,eth,eg,p ;

#if NO_MOV
	fprintf(fp,"%d %d \n",NX,NY) ;
#endif
        for(j=0;j<NY;j++) 
        for(i=0;i<NX;i++) {
                vxa = 0.5*(vx[i+1][j]+vx[i][j]) ;
                vya = 0.5*(vy[i][j+1]+vy[i][j]) ;
                vza = vz[i][j] ;
                Bxa = 0.5*(Bx[i+1][j]+Bx[i][j]) ;
                Bya = 0.5*(By[i][j+1]+By[i][j]) ;
                Bza = Bz[i][j] ;

		if(fabs(gam - 1.) > 1.e-6) p = (gam - 1.)*e[i][j] ;
		else p = cs*cs*rho[i][j] ;

#if MOVIE
		fprintf(fp,"%3d %3d %7.4g\n", i, j, rho[i][j]) ;
#endif
#if NO_MOV
		fprintf(fp,"%15.10g %15.10g ", (i-NX/2+0.5)*a*dx, 
			(j-NY/2+0.5)*a*dy) ;
		fprintf(fp,"%15.10g %15.10g ", rho[i][j],p) ;
		fprintf(fp,"%15.10g %15.10g %15.10g ", vxa,vya,vza) ;
		fprintf(fp,"%15.10g %15.10g %15.10g ", Bxa,Bya,Bza) ;
		fprintf(fp,"%15.10g\n",pot[i][j]) ;
#endif
	}
}
