
#include "decs.h"

#define MOVIE	0
#define NO_MOV	1

void dump(fp)
FILE *fp ;
{
	int i,k ;
	static int slice = NY/2 ;
	REAL vxa,vya,vza ;
	REAL Bxa,Bya,Bza,divb,p ;
	REAL eb,ek,eth,eg ;

#if NO_MOV
	fprintf(fp,"%d %d \n",NX,NZ) ;
#endif
        for(k=0;k<NZ;k++)
        for(i=0;i<NX;i++) {
                vxa = 0.5*(vx[i+1][slice][k]+vx[i][slice][k]) ;
                vya = 0.5*(vy[i][slice+1][k]+vy[i][slice][k]) ;
                vza = vz[i][slice][k] ;
                Bxa = 0.5*(bx[i+1][slice][k]+bx[i][slice][k]) ;
                Bya = 0.5*(by[i][slice+1][k]+by[i][slice][k]) ;
                Bza = bz[i][slice][k] ;

		divb = (bx[i+1][slice][k] - bx[i][slice][k])/dx
		     + (by[i][slice+1][k] - by[i][slice][k])/dy
		     + (bz[i][slice][k+1] - bz[i][slice][k])/dz ;
		
                if(fabs(gam - 1.) > 1.e-6) p = (gam - 1.)*e[i][slice][k] ;
                else p = cs*cs*rho[i][slice][k] ;
#if MOVIE
		fprintf(fp,"%3d %3d %7.4g\n", i, k, rho[i][slice][k]) ;
#endif
#if NO_MOV
		fprintf(fp,"%15.10g %15.10g ", (i-NX/2+0.5)*dx,
			(slice-NY/2+0.5)*dy) ;
		fprintf(fp,"%15.10g %15.10g ", rho[i][slice][k],p) ;
		fprintf(fp,"%15.10g %15.10g %15.10g ", vxa,vya,vza) ;
		fprintf(fp,"%15.10g %15.10g %15.10g ", Bxa,Bya,Bza) ;
		fprintf(fp,"%15.10g\n",divb) ;
#endif
	}
}
