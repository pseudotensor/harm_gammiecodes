
#include <stdio.h>
#include <math.h>

#define NX	128	/* number of zones */
#define NY	128	/* number of zones */
#define NP	8	/* number of primitive variables */

#define RHO	0	/* mnemonics for primitive vars */
#define UX	1
#define UY	2
#define UZ	3
#define BX	4
#define BY	5
#define BZ	6
#define UU	7

extern double   a_p[NX+2][NY+2][NP] ;	/* space for primitive vars */
extern double a_dqx[NX+2][NY+2][NP] ;	/* slopes */
extern double a_dqy[NX+2][NY+2][NP] ;	/* slopes */
extern double  a_Fx[NX+2][NY+2][NP] ;	/* fluxes */
extern double  a_Fy[NX+2][NY+2][NP] ;	/* fluxes */
extern double  a_ph[NX+2][NY+2][NP] ;	/* half-step primitives */

extern double (*   p)[NY+2][NP] ;
extern double (* dqx)[NY+2][NP] ;
extern double (* dqy)[NY+2][NP] ;
extern double (*  Fx)[NY+2][NP] ;
extern double (*  Fy)[NY+2][NP] ;
extern double (*  ph)[NY+2][NP] ;
extern double (* emf)[NY+2][NP] ;

extern double Lx,Ly ;

extern double cour ;
extern double cs ;
extern double gam ;
extern double dx,dy,dV ;
extern double dt ;
extern double t,tf ;

extern double DTd ;
extern double DTl ;
extern double DTi ;

/* choose boundary conditions.  Only one should be
   set to 1, all others to 0. */
#define	PER	1	/* periodic */
#define	REF	0	/* reflecting */
#define STD	0	/* set ghost zones equal to adjacent zones */

/* numerical convenience */
#define SMALL	1.e-10 

#define LOOP for(i=0;i<NX;i++)for(j=0;j<NY;j++)
#define PLOOP for(k=0;k<NP;k++)

