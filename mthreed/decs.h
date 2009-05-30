
#include <stdio.h>
#include <math.h>

#define NX	128
#define NY	128
#define NZ	32

#define SMALL	1.e-10
#define GMIN 	1.e-10
#define MIN	1.e-20	/* minimum density */

#define REAL	double

#define LINEAR   	1  /* linear theory test */
#define HDMODE1  	0  /* epicyclic mode */
#define HDMODE2  	0  /* rotationally-modified sound waves */
#define MHDMODE1 	0  /* Alfven waves */
#define MHDMODE2 	0  /* magnetosonic waves */
#define MHDMODE3 	0  /* MRI */

#define RHOCOL		0  /* circular density column */
#define STRIPE		0  /* horizontal density stripe */
#define FLOOP		0  /* B-field loop */
#define BY0		0  /* linear time dependence of By */
#define BVANLEER	1  /* use van Leer slopes for B-field shear advection */

#define SWEEPYM		1  /* mean-velocity advection */
#define LINEAR_REMAP	0  /* linear interpolation */
#define FLUX_REMAP	1  /* flux interpolation */
#define FOURIER_REMAP	0  /* fourier interpolation */

extern REAL   rhoas[NX+5][NY+5][NZ+5] ;	
extern REAL   rhotmpas [NX+5][NY+5][NZ+5] ;
extern REAL   vxas[NX+5][NY+5][NZ+5] ;
extern REAL   vyas[NX+5][NY+5][NZ+5] ;
extern REAL   vzas[NX+5][NY+5][NZ+5] ;
extern REAL   bxas[NX+5][NY+5][NZ+5] ;
extern REAL   byas[NX+5][NY+5][NZ+5] ;
extern REAL   bzas[NX+5][NY+5][NZ+5] ;
extern REAL   eas[NX+5][NY+5][NZ+5] ;
extern REAL   work1as[NX+4][NY+4][NZ+4] ;
extern REAL   work2as[NX+4][NY+4][NZ+4] ;
extern REAL   work3as[NX+4][NY+4][NZ+4] ;
extern REAL   work4as[NX+4][NY+4][NZ+4] ;
extern REAL   work5as[NX+4][NY+4][NZ+4] ;
extern REAL   work6as[NX+4][NY+4][NZ+4] ;

extern REAL (* rho)[NY+5][NZ+5] ;
extern REAL (* rhotmp)[NY+5][NZ+5] ;
extern REAL (*  vx)[NY+5][NZ+5] ;
extern REAL (*  vy)[NY+5][NZ+5] ;
extern REAL (*  vz)[NY+5][NZ+5] ;
extern REAL (*  bx)[NY+5][NZ+5] ;
extern REAL (*  by)[NY+5][NZ+5] ;
extern REAL (*  bz)[NY+5][NZ+5] ;
extern REAL (*   e)[NY+5][NZ+5] ;
extern REAL (*  work1)[NY+4][NZ+4] ;
extern REAL (*  work2)[NY+4][NZ+4] ;
extern REAL (*  work3)[NY+4][NZ+4] ;
extern REAL (*  work4)[NY+4][NZ+4] ;
extern REAL (*  work5)[NY+4][NZ+4] ;
extern REAL (*  work6)[NY+4][NZ+4] ;

extern REAL cour ;
extern REAL cs ;
extern REAL c ;
extern REAL gam ;
extern REAL dx,dy,dz ;
extern REAL dt ;
extern REAL t,tf ;
extern REAL nu_vnr ;
extern REAL nu_l ;

extern REAL Lx ;
extern REAL Ly ;
extern REAL Lz ;

extern REAL dV ;

extern REAL DTd ;
extern REAL DTi ;
extern REAL DTl ;

extern REAL drive ;
extern REAL A0,B0,B0x,B0y,B0z ;
extern REAL work,workp ;
extern REAL W,Wz,q ;

extern REAL kx,ky,kz ;
extern REAL remr[8],remi[8],norm,pert ;

extern int nstep ;

void remap_linear(REAL (*var)[NY+5][NZ+5], int i, int k, REAL shift) ;
void remap_fourier(REAL (*var)[NY+5][NZ+5], int i, int k, REAL del) ;
void remap_flux(REAL (*var)[NY+5][NZ+5], int i, int k, int shift, REAL del) ;
void consistent_transport(REAL (*var)[NY+5][NZ+5], int i, int k, REAL *massvar, REAL *specvar, REAL *fluxvar, int fr) ;
void flux_calc(REAL *var, REAL *trnsvar, REAL *mdotvar, REAL *flux, int offset) ;
void flux_update(REAL *fluxvar, REAL *flux, int offset) ;
void iremap(REAL *var, int shift) ;
int my_nint(double dum) ;
REAL slope_lim(double y1,double y2,double y3) ;

#define LOOP for(i=0;i<NX;i++)for(j=0;j<NY;j++)for(k=0;k<NZ;k++)
#define LOOPP(a,b,c,d,e,f) for(i=0-(a);i<NX+(b);i++)for(j=0-(c);j<NY+(d);j++)for(k=0-(e);k<NZ+(f);k++)
/* find largest axis */
#define NN     ( (NX >= NY && NX >= NZ) ? NX : ( (NY >= NX && NY >= NZ) ? NY : NZ ) )

#define ANTISYMMETRIC	-1
#define SYMMETRIC	1

