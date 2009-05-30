
#include "decs.h"

/** more Physics **/

void Utoprim(double *U, double *pr)
{
	double u1,u2,u3,u4,u5,u6,u7,u8 ;

	u1 = U[0] ;
	u2 = U[1] ;
	u3 = U[2] ;
	u4 = U[3] ;
	u5 = U[4] ;
	u6 = U[5] ;
	u7 = U[6] ;
	u8 = U[7] ;

	pr[RHO] = u1 ;
	pr[UX] = u2/u1 ;
	pr[UY] = u3/u1 ;
	pr[UZ] = u4/u1 ;
	pr[BX] = u5 ;
	pr[BY] = u6 ;
	pr[BZ] = u7 ;
	pr[UU] = u8 - 0.5*(u5*u5 + u6*u6 + u7*u7 + (u2*u2 + u3*u3 + u4*u4)/u1) ;

}

void primtoxflux(double *pr, double *f)
{
	double r,vx,vy,vz,bx,by,bz,u ;

	r = pr[RHO] ;
	vx = pr[UX] ;
	vy = pr[UY] ;
	vz = pr[UZ] ;
	bx = pr[BX] ;
	by = pr[BY] ;
	bz = pr[BZ] ;
	u = pr[UU] ;

	f[0] = r*vx ;
	f[1] = r*vx*vx + (gam - 1)*u + 0.5*(by*by + bz*bz - bx*bx) ;
	f[2] = -bx*by + r*vx*vy ;
	f[3] = -bx*bz + r*vx*vz ;
	f[4] = 0. ;
	f[5] = by*vx - bx*vy ;
	f[6] = bz*vx - bx*vz ;
	f[7] = vx*(r*0.5*(vx*vx + vy*vy + vz*vz) + by*by + bz*bz + gam*u) 
		- bx*(by*vy + bz*vz) ;
}

void primtoyflux(double *pr, double *f)
{
	double r,vx,vy,vz,bx,by,bz,u ;

	r = pr[RHO] ;
	vx = pr[UX] ;
	vy = pr[UY] ;
	vz = pr[UZ] ;
	bx = pr[BX] ;
	by = pr[BY] ;
	bz = pr[BZ] ;
	u = pr[UU] ;

	f[0] = r*vy ;
	f[1] = -by*bx + r*vy*vx ;
	f[2] = r*vy*vy + (gam - 1)*u + 0.5*(-by*by + bz*bz + bx*bx) ;
	f[3] = -by*bz + r*vy*vz ;
	f[4] = bx*vy - by*vx ;
	f[5] = 0. ;
	f[6] = bz*vy - by*vz ;
	f[7] = vy*(r*0.5*(vx*vx + vy*vy + vz*vz) + bx*bx + bz*bz + gam*u) 
		- by*(bx*vx + bz*vz) ;

}

void primtoU(double *pr, double *U)
{
	double r,vx,vy,vz,bx,by,bz,u ;

	r = pr[RHO] ;
	vx = pr[UX] ;
	vy = pr[UY] ;
	vz = pr[UZ] ;
	bx = pr[BX] ;
	by = pr[BY] ;
	bz = pr[BZ] ;
	u = pr[UU] ;

	U[0] = r ;
	U[1] = r*vx ;
	U[2] = r*vy ;
	U[3] = r*vz ;
	U[4] = bx ;
	U[5] = by ;
	U[6] = bz ;
	U[7] = 0.5*(bx*bx + by*by + bz*bz) + u + 0.5*r*(vx*vx + vy*vy + vz*vz) ;
}

