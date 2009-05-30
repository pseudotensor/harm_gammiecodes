
/* shock tube test */

#include "decs.h"

void init()
{
	int i,j ;
	void aupdate() ;
	void bound_var() ;
	void init_spec() ;
	double ranc() ;
	void set_arrays() ;
	void zero_arrays() ;
	double X,Y,R ;
	double w,K,XP,YP,amp,sig2 ;
	double beta,nj ;
	int seed ;

	seed = 1 ;
	ampf = 0.0 ;
	beta = 0. ;
	nj = 0. ;
	a = 1. ;

	/* fire up random number generator */
	ranc(seed) ;

	/* assume throughout that G = rho = L = 1 */
	cour = 0.5 ;
	cs = 1. ;
	gam = 5./3. ;
	Lx = 2.*M_PI ;
	Ly = 2.*M_PI ;
	dx = Lx/NX ;
	dy = Ly/NY ;
	dt = 0.1 ;	/* needs initial value */
	t = 0. ;
	tf = 4.0 ;
	nu_vnr = 1.0 ;
	nu_l = 0.0;

	/*
	G = M_PI*nj*nj ;
	*/
	G = 0. ;

	tdrive = 1.e9 ;
	dtdrive = 1.e9 ;
	kpk = 4. ;

	DTd = 0.1 ;	/* dumping frequency */
	DTi = 0.05 ;	/* image frequency */
	DTl = 0.01 ;	/* dumping frequency */
	res = 0. ;
	nu_sh = 0. ;

	tdecay = 1.e15 ;

	/* set up arrays,details */
	set_arrays() ;
	aupdate() ;

	if(ampf < 0.)
		init_spec() ;
	else
		zero_arrays() ;

	amp = 0.04 ;
	fprintf(stderr,"%g %g %g\n",cs,amp) ;
        for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		X = (i + .5)*dx - .5*Lx ;
		XP = i*dx - .5*Lx ;
		Y = (j + .5)*dy - .5*Ly ;
		YP = j*dy - .5*Ly ;

		Bx[i][j] = -sin(Y + M_PI) ;
		By[i][j] = sin(2.*(X + M_PI)) ;
		Bz[i][j] = 0. ;

		vx[i][j] = -sin(Y + M_PI) ;
		vy[i][j] = sin(X + M_PI) ;
		vz[i][j] = 0. ;

		rho[i][j] = 25./9. ;

		/*
		e[i][j] += cs*cs ;
		*/
		e[i][j] = (5./3.)/(gam - 1.) ;

		pot[i][j] = 0. ;
	}
	/* enforce boundary conditions */
	bound_var(Bx) ;
	bound_var(By) ;
	bound_var(Bz) ;
	bound_var(vx) ;
	bound_var(vy) ;
	bound_var(vz) ;
	bound_var(rho) ;
	bound_var(e) ;
	bound_var(pot) ;
}
