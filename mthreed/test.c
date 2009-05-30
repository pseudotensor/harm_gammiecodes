
#include <math.h>
#include <stdio.h>

float vx[6][6][6],rho[6][6][6] ;

int main()
{
	int i,j,k ;
	float Z ;

	Z = 1. ;

	for(i=1;i<=4;i++)
	for(j=1;j<=4;j++) 
	for(k=1;k<=4;k++) {
		fprintf(stdout,"%d %d %d\n",i,j,k) ;

		vx[i][j][k] = 1. ;

		fprintf(stdout,"%d %d %d\n",i,j,k) ;
		
		rho[i][j][k] = exp(-0.5*Z*Z) ;

		fprintf(stdout,"%d %d %d\n",i,j,k) ;
	}

	return(0) ;
}

