
/* initialize all variables */

#include "decs.h"

void init()
{
	int i,j,k ;
	void bound_var(REAL (* var)[NY+5][NZ+5], REAL tcurr, int sym) ;
	double ranc() ;
	void array_init() ;
	REAL X,Y,Z,R ;
	REAL w,K,XP,YP,ZP,amp,sig2 ;
	REAL b,discr,omega,omega2,kappa2,va,k0 ;
	REAL (*Az)[NY+4][NZ+4] ;
	REAL eps ;

	/* here are the real adjustable parameters */
	Lx = 16. ;
	Ly = 16*M_PI ;
	Lz = 2. ;
	A0 = 0.5 ;
	drive = 0. ;
	W = 1. ;
	Wz = 0. ;

	/* here are the seldom adjusted parameters */
	q = 1.5 ;
	kappa2 = 2.*(2. - q)*W*W ;
	cs = 1.0 ;
	c = 10.*cs ;
	gam = 1. ;

	/* numerical parameters */
	cour = 0.5 ;
	nu_vnr = 1. ;
	nu_l = 0. ;

	/* some initializations */
	dx = Lx/NX ;
	dy = Ly/NY ;
	dz = Lz/NZ ;
	dV = dx*dy*dz ;
	dt = 0.01 ;  /* needs initial value */
	t = 0. ;

	/* zero arrays & shift pointers */
	array_init() ;

	/* problem set up */
	B0z = sqrt(15.)/(32.*M_PI) ;
	tf = 1000. ;
	eps = 0.01 ;

	/* set primitive variables */
	LOOP {
		X = (i + 0.5)*dx - 0.5*Lx ;
		XP = i*dx - 0.5*Lx ;
		Y = (j + 0.5)*dy - 0.5*Ly ;
		YP = j*dy - 0.5*Ly ;
		Z = (k + 0.5)*dz - 0.5*Lz ;
		ZP = k*dz - 0.5*Lz ;
		
		vx[i][j][k] = 0. + eps*ranc(0) ;
		vy[i][j][k] = 0. + eps*ranc(0);
		vz[i][j][k] = 0. + eps*ranc(0);
		
		bx[i][j][k] = 0. ;
		by[i][j][k] = B0z ;
		bz[i][j][k] = 0. ;
		
		rho[i][j][k] = 1. ;
	}

	DTd = tf/20. ;	/* dumping frequency */
	DTi = tf/1000. ; /* image frequency */
	DTl = tf/4000. ; /* log frequency */

	/* enforce boundary conditions */
	bound_var(rho,t,SYMMETRIC) ;
	bound_var(vx,t,SYMMETRIC) ;
	bound_var(vy,t,SYMMETRIC) ;
	bound_var(vz,t,ANTISYMMETRIC) ;

	bound_var(by,t,SYMMETRIC) ;
	bound_var(bz,t,ANTISYMMETRIC) ;
        bound_var(bx,t,SYMMETRIC) ;

	bound_var(e,t,SYMMETRIC) ;
}
