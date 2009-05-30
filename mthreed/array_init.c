
#include "decs.h"

void array_init()
{
	int i,j,k ;

        /* pointer shifting */
        rho =  (REAL (*) [NY+5][NZ+5])(& ( rhoas[3][3][3])) ;
#if FLUX_REMAP
	rhotmp =  (REAL (*) [NY+5][NZ+5])(& ( rhotmpas[3][3][3])) ;
#endif
        vx =   (REAL (*) [NY+5][NZ+5])(& (  vxas[3][3][3])) ;
        vy =   (REAL (*) [NY+5][NZ+5])(& (  vyas[3][3][3])) ;
        vz =   (REAL (*) [NY+5][NZ+5])(& (  vzas[3][3][3])) ;
        bx =   (REAL (*) [NY+5][NZ+5])(& (  bxas[3][3][3])) ;
        by =   (REAL (*) [NY+5][NZ+5])(& (  byas[3][3][3])) ;
        bz =   (REAL (*) [NY+5][NZ+5])(& (  bzas[3][3][3])) ;
	e  =   (REAL (*) [NY+5][NZ+5])(& (   eas[3][3][3])) ;
        work1 =   (REAL (*) [NY+4][NZ+4])(& (  work1as[2][2][2])) ;
        work2 =   (REAL (*) [NY+4][NZ+4])(& (  work2as[2][2][2])) ;
        work3 =   (REAL (*) [NY+4][NZ+4])(& (  work3as[2][2][2])) ;
        work4 =   (REAL (*) [NY+4][NZ+4])(& (  work4as[2][2][2])) ;
        work5 =   (REAL (*) [NY+4][NZ+4])(& (  work5as[2][2][2])) ;
        work6 =   (REAL (*) [NY+4][NZ+4])(& (  work6as[2][2][2])) ;

	/* zero arrays */
	LOOPP(3,2,3,2,3,2)  rho[i][j][k] = 0. ;

	LOOPP(2,2,2,2,2,2) {
		 vx[i][j][k] = 0. ;
		 vy[i][j][k] = 0. ;
		 vz[i][j][k] = 0. ;
		 bx[i][j][k] = 0. ;
		 by[i][j][k] = 0. ;
		 bz[i][j][k] = 0. ;
		 e[i][j][k] = 0. ;

		 work1[i][j][k] = 0. ;
		 work2[i][j][k] = 0. ;
		 work3[i][j][k] = 0. ;
		 work4[i][j][k] = 0. ;
		 work5[i][j][k] = 0. ;
		 work6[i][j][k] = 0. ;
	}
	
#if FLUX_REMAP
	LOOPP(3,2,2,2,2,2) rhotmp[i][j][k] = 0. ; 
#endif
}
