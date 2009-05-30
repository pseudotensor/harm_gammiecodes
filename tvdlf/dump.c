
#include "decs.h"

void dump(fp)
FILE *fp ;
{
	int i,j,k ;

	fprintf(fp,"%d %d %10.5g %10.5g %10.5g\n",NX,NY,Lx,Ly,t) ;

	LOOP {	
		fprintf(fp,"%15.7g %15.7g",(i-NX/2+0.5)*dx,(j-NY/2+0.5)*dy) ;
		PLOOP fprintf(fp,"%15.7g ",p[i][j][k]) ;
		fprintf(fp,"\n") ;
	}

}
