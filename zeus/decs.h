
#include <stdio.h>
#include <math.h>

extern double **rho ;	
extern double **e ;
extern double **vx ;
extern double **vy ;
extern double **vz ;
extern double **Bx ;
extern double **By ;
extern double **Bz ;
extern double **pot ;
extern double **fl ;
extern double **work1 ;
extern double **work2 ;
extern double **work3 ;
extern double **work4 ;
extern double **work5 ;
extern double **work6 ;
extern double **work7 ;
extern double **work8 ;
extern double **work9 ;
extern double **work10 ;

extern double cour ;
extern double cs ;
extern double gam ;
extern double dx,dy ;
extern double dt ;
extern double t,tf ;
extern double nu_vnr ;
extern double nu_l ;
extern double G ;

extern int NX ;
extern int NY ;

extern double Lx ;
extern double Ly ;

extern double DTd ;
extern double DTi ;
extern double DTl ;
extern double res ;
extern double nu_sh ;

extern double a ;
extern double da ;
extern double tdecay ;
extern double ampf ;
extern double tdrive ;
extern double dtdrive ;
extern double kpk ;

#define SMALL	1.e-10
#define GMIN 	1.e-10
#define MIN	1.e-20	/* minimum density */

#define DONOR	0	/* use donor cell or van leer? */
#define VANLEER	1

