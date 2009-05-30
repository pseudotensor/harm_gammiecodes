
#include "decs.h"

void step_visc()
{
	/* could minimize storage by only using one vis? variable
	   and looping over direction. */
	static REAL visxa[NX+4][NY+4][NZ+4] ;
	static REAL visya[NX+4][NY+4][NZ+4] ;
	static REAL visza[NX+4][NY+4][NZ+4] ;
	static REAL (*visx)[NY+4][NZ+4] ;
	static REAL (*visy)[NY+4][NZ+4] ;
	static REAL (*visz)[NY+4][NZ+4] ;
	REAL dvx,dvy,dvz,qlx,qvnrx,qly,qvnry,qlz,qvnrz ;
	REAL rhoa,rhoi ;
	int i,j,k ;
	void bound_var(REAL (*var)[NY+5][NZ+5], REAL tcurr, int sym) ;

	/* pointer shifting */
	visx = (REAL (*) [NY+4][NZ+4])( & (visxa[2][2][2])) ;
	visy = (REAL (*) [NY+4][NZ+4])( & (visya[2][2][2])) ;
	visz = (REAL (*) [NY+4][NZ+4])( & (visza[2][2][2])) ;

	/* find viscous stresses */
	LOOPP(1,0,1,0,1,0) {
		/* del v, at zone center */
		dvx = vx[i+1][j][k] - vx[i][j][k] ;
		dvy = vy[i][j+1][k] - vy[i][j][k] ;
		dvz = vz[i][j][k+1] - vz[i][j][k] ;

		/* linear viscosity */
		if(fabs(gam - 1.) > 1.e-6) cs = sqrt(gam*(gam-1.)*e[i][j][k]/rho[i][j][k]) ;
		rhoi = rho[i][j][k] ;
		qlx = -nu_l*rhoi*cs*dvx ;
		qly = -nu_l*rhoi*cs*dvy ;
		qlz = -nu_l*rhoi*cs*dvz ;

		/* von neumann,richtmyer viscosity */
		qvnrx = (dvx < 0.) ? nu_vnr*rhoi*dvx*dvx : 0. ;
		qvnry = (dvy < 0.) ? nu_vnr*rhoi*dvy*dvy : 0. ;
		qvnrz = (dvz < 0.) ? nu_vnr*rhoi*dvz*dvz : 0. ;

		visx[i][j][k] = qlx + qvnrx ;
		visy[i][j][k] = qly + qvnry ;
		visz[i][j][k] = qlz + qvnrz ;
	}
	/* don't bound! */

	/* update velocity, internal energy */
	if(fabs(gam - 1.) > 1.e-6) {
	LOOP {
		e[i][j][k] +=
			   - dt*visx[i][j][k]*(vx[i+1][j][k] - vx[i][j][k])/dx
			   - dt*visy[i][j][k]*(vy[i][j+1][k] - vy[i][j][k])/dy
			   - dt*visz[i][j][k]*(vz[i][j][k+1] - vz[i][j][k])/dz ;
	}
	bound_var(e,t,SYMMETRIC) ;
	}

	LOOP {
		rhoi = rho[i][j][k] ;

		rhoa = 0.5*(rhoi+rho[i-1][j][k]) ;
		vx[i][j][k] += -dt*(visx[i][j][k]-visx[i-1][j][k])/(dx*rhoa) ;

		rhoa = 0.5*(rhoi+rho[i][j-1][k]) ;
		vy[i][j][k] += -dt*(visy[i][j][k]-visy[i][j-1][k])/(dy*rhoa) ;

		rhoa = 0.5*(rhoi+rho[i][j][k-1]) ;
		vz[i][j][k] += -dt*(visz[i][j][k]-visz[i][j][k-1])/(dz*rhoa) ;
	}
	bound_var(vx,t,SYMMETRIC) ;
	bound_var(vy,t,SYMMETRIC) ;
	bound_var(vz,t,ANTISYMMETRIC) ;
}
