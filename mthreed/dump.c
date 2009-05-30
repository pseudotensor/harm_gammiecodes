
#include "decs.h"

void dump(fp)
FILE *fp ;
{
	int i,j,k ;
	REAL vxa,vya,vza ;
	REAL Bxa,Bya,Bza,divb,p ;
	REAL eb,ek,eth,eg ;

	/* first line in dump file */
	fprintf(fp,"%d %d %d %g %g %g %g\n",NX,NY,NZ,t,Lx,Ly,Lz) ;

        for(i=0;i<NX;i++) 
        for(j=0;j<NY;j++) 
        for(k=0;k<NZ;k++) 
	{
		/* get everything centered on the zone */
                vxa = 0.5*(vx[i+1][j][k]+vx[i][j][k]) ;
                vya = 0.5*(vy[i][j+1][k]+vy[i][j][k]) ;
                vza = 0.5*(vz[i][j][k+1]+vz[i][j][k]) ;

                Bxa = 0.5*(bx[i+1][j][k]+bx[i][j][k]) ;
                Bya = 0.5*(by[i][j+1][k]+by[i][j][k]) ;
                Bza = 0.5*(bz[i][j][k+1]+bz[i][j][k]) ;

		fprintf(fp,"%d %d %d", i, j, k) ;
		fprintf(fp,"%15.10g ", rho[i][j][k]) ;
		fprintf(fp,"%15.10g %15.10g %15.10g ", vxa,vya,vza) ;
		fprintf(fp,"%15.10g %15.10g %15.10g ", Bxa,Bya,Bza) ;
	}

}
