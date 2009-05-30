
#include "decs.h"

/* pressure and gravity step */
void step_pg()
{
	REAL rhoa,XP,ZP,vxave,vyave ;
	REAL bxa,bxp,bxm,bya,byp,bym,bza,bzp,bzm ;
	int i,j,k ;
	void bound_var(REAL (*var)[NY+5][NZ+5], REAL tcurr, int sym) ;

	/* gas pressure gradient acceleration */
	if(fabs(gam - 1.) > 1.e-6) {
	LOOP {
		rhoa = 0.5*(rho[i][j][k]+rho[i-1][j][k]) ;
		vx[i][j][k] += -(dt/(2.*dx*rhoa))*(gam-1.)*(e[i][j][k]-e[i-1][j][k]) ;

		rhoa = 0.5*(rho[i][j][k]+rho[i][j][k-1]) ;
		vz[i][j][k] += -(dt/(2.*dz*rhoa))*(gam-1.)*(e[i][j][k]-e[i][j][k-1]) ;
	}
	}
	else {
	LOOP {
		rhoa = 0.5*(rho[i][j][k]+rho[i-1][j][k]) ;
		vx[i][j][k] += -(dt/(2.*dx*rhoa))*cs*cs*(rho[i][j][k]-rho[i-1][j][k]) ;

		rhoa = 0.5*(rho[i][j][k]+rho[i][j][k-1]) ;
		vz[i][j][k] += -(dt/(2.*dz*rhoa))*cs*cs*(rho[i][j][k]-rho[i][j][k-1]) ;
	}
	}

	/* magnetic pressure gradient acceleration */
	LOOP {
		rhoa = 0.5*(rho[i][j][k]+rho[i-1][j][k]) ;
		byp = 0.5*(by[i][j][k] + by[i][j+1][k]) ;
		bym = 0.5*(by[i-1][j][k] + by[i-1][j+1][k]) ;
		bya = 0.5*(byp + bym) ;

		bzp = 0.5*(bz[i][j][k] + bz[i][j][k+1]) ;
		bzm = 0.5*(bz[i-1][j][k] + bz[i-1][j][k+1]) ;
		bza = 0.5*(bzp + bzm) ;
		vx[i][j][k] += -(dt/(2.*dx*rhoa))*(bya*(byp - bym) + bza*(bzp - bzm)) ;

		rhoa = 0.5*(rho[i][j][k] + rho[i][j][k-1]) ;
		bxp = 0.5*(bx[i][j][k] + bx[i+1][j][k]) ;
		bxm = 0.5*(bx[i][j][k-1] + bx[i+1][j][k-1]) ;
		bxa = 0.5*(bxp + bxm) ;

		byp = 0.5*(by[i][j][k] + by[i][j+1][k]) ;
		bym = 0.5*(by[i][j][k-1] + by[i][j+1][k-1]) ;
		bya = 0.5*(byp + bym) ;
		vz[i][j][k] += -(dt/(2.*dz*rhoa))*(bxa*(bxp - bxm) + bya*(byp - bym)) ;
	}

        /* vertical tidal force */
	LOOP {
                ZP = k*dz - 0.5*Lz ;

                vz[i][j][k] += (dt/2.)*(-Wz*Wz*ZP) ;
        }

        /* Coriolis force, tidal acc */
	if(W!=0.) LOOP {
                vyave = 0.25*(vy[i][j][k] + vy[i-1][j][k] +
			vy[i][j+1][k] + vy[i-1][j+1][k]) ;
                vx[i][j][k] += (dt/2.)*(2.*W*vyave) ;
#if SWEEPYM
#else
		XP = i*dx - 0.5*Lx ;
		vx[i][j][k] += (dt/2.)*2.*q*W*W*XP ;
#endif
        }
	bound_var(vx,t,SYMMETRIC) ;

      	/* gas pressure gradient acceleration */
	if(fabs(gam - 1.) > 1.e-6) {
	LOOP {
		rhoa = 0.5*(rho[i][j][k]+rho[i][j-1][k]) ;
		vy[i][j][k] += -(dt/(dy*rhoa))*(gam-1.)*(e[i][j][k]-e[i][j-1][k]) ;
	}
	}
	else {
	LOOP {
		rhoa = 0.5*(rho[i][j][k]+rho[i][j-1][k]) ;
		vy[i][j][k] += -(dt/(dy*rhoa))*cs*cs*(rho[i][j][k]-rho[i][j-1][k]) ;
	}
	}

	/* magnetic pressure gradient acceleration */
	LOOP {
		rhoa = 0.5*(rho[i][j][k] + rho[i][j-1][k]) ;
		bxp = 0.5*(bx[i][j][k] + bx[i+1][j][k]) ;
		bxm = 0.5*(bx[i][j-1][k] + bx[i+1][j-1][k]) ;
		bxa = 0.5*(bxp + bxm) ;

		bzp = 0.5*(bz[i][j][k] + bz[i][j][k+1]) ;
		bzm = 0.5*(bz[i][j-1][k] + bz[i][j-1][k+1]) ;
		bza = 0.5*(bzp + bzm) ;
		vy[i][j][k] += -(dt/(dy*rhoa))*(bxa*(bxp - bxm) + bza*(bzp - bzm)) ;
	}

	if(W!=0.) LOOP {
                vxave = 0.25*(vx[i][j][k] + vx[i+1][j][k] +
			vx[i][j-1][k] + vx[i+1][j-1][k]) ;
                vy[i][j][k] += dt*(-2.*W*vxave) ;
#if SWEEPYM
		vy[i][j][k] += dt*q*W*vxave ;
#endif
        }
        bound_var(vy,t,SYMMETRIC) ;
 
	/* gas pressure gradient acceleration */
	if(fabs(gam - 1.) > 1.e-6) {
	LOOP {
		rhoa = 0.5*(rho[i][j][k]+rho[i-1][j][k]) ;
		vx[i][j][k] += -(dt/(2.*dx*rhoa))*(gam-1.)*(e[i][j][k]-e[i-1][j][k]) ;

		rhoa = 0.5*(rho[i][j][k]+rho[i][j][k-1]) ;
		vz[i][j][k] += -(dt/(2.*dz*rhoa))*(gam-1.)*(e[i][j][k]-e[i][j][k-1]) ;
	}
	}
	else {
	LOOP {
		rhoa = 0.5*(rho[i][j][k]+rho[i-1][j][k]) ;
		vx[i][j][k] += -(dt/(2.*dx*rhoa))*cs*cs*(rho[i][j][k]-rho[i-1][j][k]) ;

		rhoa = 0.5*(rho[i][j][k]+rho[i][j][k-1]) ;
		vz[i][j][k] += -(dt/(2.*dz*rhoa))*cs*cs*(rho[i][j][k]-rho[i][j][k-1]) ;
	}
	}

	/* magnetic pressure gradient acceleration */
	LOOP {
		rhoa = 0.5*(rho[i][j][k]+rho[i-1][j][k]) ;
		byp = 0.5*(by[i][j][k] + by[i][j+1][k]) ;
		bym = 0.5*(by[i-1][j][k] + by[i-1][j+1][k]) ;
		bya = 0.5*(byp + bym) ;

		bzp = 0.5*(bz[i][j][k] + bz[i][j][k+1]) ;
		bzm = 0.5*(bz[i-1][j][k] + bz[i-1][j][k+1]) ;
		bza = 0.5*(bzp + bzm) ;
		vx[i][j][k] += -(dt/(2.*dx*rhoa))*(bya*(byp - bym) + bza*(bzp - bzm)) ;

		rhoa = 0.5*(rho[i][j][k] + rho[i][j][k-1]) ;
		bxp = 0.5*(bx[i][j][k] + bx[i+1][j][k]) ;
		bxm = 0.5*(bx[i][j][k-1] + bx[i+1][j][k-1]) ;
		bxa = 0.5*(bxp + bxm) ;

		byp = 0.5*(by[i][j][k] + by[i][j+1][k]) ;
		bym = 0.5*(by[i][j][k-1] + by[i][j+1][k-1]) ;
		bya = 0.5*(byp + bym) ;
		vz[i][j][k] += -(dt/(2.*dz*rhoa))*(bxa*(bxp - bxm) + bya*(byp - bym)) ;
	}

        /* vertical tidal force */
	LOOP {
                ZP = k*dz - 0.5*Lz ;

                vz[i][j][k] += (dt/2.)*(-Wz*Wz*ZP) ;
        }
        bound_var(vz,t,ANTISYMMETRIC) ;

        /* Coriolis force, tidal acc */
	if(W!=0.) LOOP {
                vyave = 0.25*(vy[i][j][k] + vy[i-1][j][k] +
			vy[i][j+1][k] + vy[i-1][j+1][k]) ;
                vx[i][j][k] += (dt/2.)*(2.*W*vyave) ;
#if SWEEPYM
#else
		XP = i*dx - 0.5*Lx ;
		vx[i][j][k] += (dt/2.)*2.*q*W*W*XP ;
#endif
        }
        bound_var(vx,t,SYMMETRIC) ;
}
