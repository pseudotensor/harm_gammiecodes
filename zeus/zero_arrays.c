
#include "decs.h"

void zero_arrays()
{
	int i,j ;

	for(i=-1;i<=NX;i++)
	for(j=-1;j<=NY;j++) {
		rho[i][j] = 0. ;
		e[i][j] = 0. ;
		vx[i][j] = 0. ;
		vy[i][j] = 0. ;
		vz[i][j] = 0. ;
		Bx[i][j] = 0. ;
		By[i][j] = 0. ;
		Bz[i][j] = 0. ;
		pot[i][j] = 0. ;
		fl[i][j] = 0. ;
		work1[i][j] = 0. ;
		work2[i][j] = 0. ;
		work3[i][j] = 0. ;
		work4[i][j] = 0. ;
		work5[i][j] = 0. ;
		work6[i][j] = 0. ;
		work7[i][j] = 0. ;
		work8[i][j] = 0. ;
		work9[i][j] = 0. ;
		work10[i][j] = 0. ;
	}
}
