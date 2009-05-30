
#include "decs.h"

/* Toth/Rusanov/Yee  method */
#define FAC	0.5	

void step_ch()
{
	int i,j,k ;
	double Dqm,Dqp,Dqc,s ;
	double (* Fx_ct)[NY+2][NP] ;
	double (* Fy_ct)[NY+2][NP] ;
	double cmax,cs2x,cs2y,va2x,va2y,cmsx,cmsy,vx,vy ;
	static double cx[NX][NY],cy[NX][NY] ;
	static double c = 1. ;
	static double Uh[NP],U[NP] ;
	static double pix[NP],plx[NP],prx[NP],Flx[NP],Frx[NP],Ulx[NP],Urx[NP] ;
	static double piy[NP],ply[NP],pry[NP],Fly[NP],Fry[NP],Uly[NP],Ury[NP] ;
	void bounds() ;
	void primtoU(double *pb, double *Ub) ;
	void Utoprim(double *Ua, double *pa) ;
	void primtoxflux(double *pa, double *Fa) ;
	void primtoyflux(double *pa, double *Fa) ;

	/** set timestep **/
	dt = cour*dx/1. ;

        /* don't step beyond end of run */
        if(t + dt > tf) dt = tf - t ;

	bounds() ;
#if 0
	/* evaluate Woodward slopes of primitive variables */
	LOOP {
		PLOOP {
			Dqm = 2.0*(p[i][j][k] - p[i-1][j][k]) ;
			Dqp = 2.0*(p[i+1][j][k] - p[i][j][k]) ;
			Dqc = 0.5*(p[i+1][j][k] - p[i-1][j][k]) ;
			s = Dqm*Dqp ;
			if(s <= 0.) dqx[i][j][k] = 0. ;
			else {
				if(fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
					dqx[i][j][k] = Dqm ;
				else if(fabs(Dqp) < fabs(Dqc))
					dqx[i][j][k] = Dqp ;
				else
					dqx[i][j][k] = Dqc ;
			}
		}
		PLOOP {
			Dqm = 2.0*(p[i][j][k] - p[i][j-1][k]) ;
			Dqp = 2.0*(p[i][j+1][k] - p[i][j][k]) ;
			Dqc = 0.5*(p[i][j+1][k] - p[i][j-1][k]) ;
			s = Dqm*Dqp ;
			if(s <= 0.) dqy[i][j][k] = 0. ;
			else {
				if(fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
					dqy[i][j][k] = Dqm ;
				else if(fabs(Dqp) < fabs(Dqc))
					dqy[i][j][k] = Dqp ;
				else
					dqy[i][j][k] = Dqc ;
			}
		}
	}
#endif
#if 1
	/* evaluate vanleer slopes */
	LOOP {
		PLOOP {
			Dqm = p[i][j][k] - p[i-1][j][k] ;
			Dqp = p[i+1][j][k] - p[i][j][k] ;
			s = Dqm*Dqp ;
			if(s <= 0.) dqx[i][j][k] = 0. ;
			else dqx[i][j][k] = 2.*s/(Dqm+Dqp) ;
		}
		PLOOP {
			Dqm = p[i][j][k] - p[i][j-1][k] ;
			Dqp = p[i][j+1][k] - p[i][j][k] ;
			s = Dqm*Dqp ;
			if(s <= 0.) dqy[i][j][k] = 0. ;
			else dqy[i][j][k] = 2.*s/(Dqm+Dqp) ;
		}
	}
#endif
	bounds() ;

	/* evaluate timestep */
        cmax = 0. ;
        LOOP {
                PLOOP {
                        pix[k] = 0.5*(
				p[i-1][j][k] + 0.5*dqx[i-1][j][k] +
                                p[i][j][k]   - 0.5*dqx[i][j][k]) ;
                        piy[k] = 0.5*(
				p[i][j-1][k] + 0.5*dqy[i][j-1][k] +
                                p[i][j][k]   - 0.5*dqy[i][j][k]) ;
                }
		// MARK
		if(pix[7]<0.){
		  fprintf(stderr,"%d %d pix[7]: %15.10g\n",i,j,pix[7]);
		}
		if(piy[7]<0.){
		  fprintf(stderr,"%d %d piy[7]: %15.10g\n",i,j,piy[7]);
		}
		if(pix[0]<0.){
		  fprintf(stderr,"%d %d pix[0]: %15.10g\n",i,j,pix[0]);
		}
		if(piy[0]<0.){
		  fprintf(stderr,"%d %d piy[0]: %15.10g\n",i,j,piy[0]);
		}
                /* evaluate speeds */
                cs2x = gam*(gam - 1.)*pix[UU]/pix[RHO] ;
                cs2y = gam*(gam - 1.)*piy[UU]/piy[RHO] ;
                va2x = (pix[BX]*pix[BX] + pix[BY]*pix[BY] + pix[BZ]*pix[BZ])/
				pix[RHO] ;
                va2y = (piy[BX]*piy[BX] + piy[BY]*piy[BY] + piy[BZ]*piy[BZ])/
				piy[RHO] ;
                cmsx = sqrt(cs2x + va2x) ;
                cmsy = sqrt(cs2y + va2y) ;
		vx = fabs(pix[UX]) ;
		vy = fabs(piy[UY]) ;

		cx[i][j] = vx + cmsx ;
		cy[i][j] = vy + cmsy ;

                if(cmsx > cmax) cmax = cmsx ;
                if(cmsy > cmax) cmax = cmsy ;
        }
        dt = cour*dx/cmax ;
        /* don't step beyond end of run */
        if(t + dt > tf) dt = tf - t ;


	/* perform Hancock predictor half-step */
	LOOP {
		PLOOP {
			plx[k] = p[i][j][k] - 0.5*dqx[i][j][k] ;
			prx[k] = p[i][j][k] + 0.5*dqx[i][j][k] ;
			ply[k] = p[i][j][k] - 0.5*dqy[i][j][k] ;
			pry[k] = p[i][j][k] + 0.5*dqy[i][j][k] ;
		}

		primtoxflux(plx,Flx) ;
		primtoxflux(prx,Frx) ;

		primtoyflux(ply,Fly) ;
		primtoyflux(pry,Fry) ;

		primtoU(p[i][j],Uh) ;

		PLOOP {
			Uh[k] -= 0.5*dt*( (Frx[k] - Flx[k])/dx 
				+ (Fry[k] - Fly[k])/dy ) ;
			ph[i][j][k] = p[i][j][k] ;
		}

		Utoprim(Uh,ph[i][j]) ;

	}
	bounds() ;

	/* evaluate zone boundary fluxes at half-step */
	LOOP {
		PLOOP {
			pix[k] = 0.5*(ph[i-1][j][k] + 0.5*dqx[i-1][j][k] +
			   	      ph[i][j][k]   - 0.5*dqx[i][j][k]) ;
			piy[k] = 0.5*(ph[i][j-1][k] + 0.5*dqy[i][j-1][k] +
			   	      ph[i][j][k]   - 0.5*dqy[i][j][k]) ;
		}
		primtoxflux(pix, Fx[i][j]) ;
		primtoyflux(piy, Fy[i][j]) ;
	}

	/* evaluate diffusive correction to flux */
	LOOP {
		PLOOP {
			plx[k] = ph[i-1][j][k] + 0.5*dqx[i-1][j][k] ;
			prx[k] = ph[i][j][k]   - 0.5*dqx[i][j][k] ;

			ply[k] = ph[i][j-1][k] + 0.5*dqy[i][j-1][k] ;
			pry[k] = ph[i][j][k]   - 0.5*dqy[i][j][k] ;
		}
		primtoU(plx,Ulx) ;
		primtoU(prx,Urx) ;
		primtoU(ply,Uly) ;
		primtoU(pry,Ury) ;
		
		PLOOP {
			Fx[i][j][k] += FAC*cx[i][j]*(Ulx[k] - Urx[k]) ;
			Fy[i][j][k] += FAC*cy[i][j]*(Uly[k] - Ury[k]) ;
		}
	}
	bounds() ;

#if 1
	/* replace B-field fluxes with flux-ct symmetrized fluxes */
	Fx_ct = dqx ;
	Fy_ct = dqy ;
	LOOP {

		Fx_ct[i][j][BX] = 0. ;
		Fy_ct[i][j][BX] = 0.125*( 
			2.*Fy[i][j][BX]
			+ Fy[i+1][j][BX]
			+ Fy[i-1][j][BX]
			- Fx[i][j][BY]
			- Fx[i][j-1][BY]
			- Fx[i+1][j][BY]
			- Fx[i+1][j-1][BY] ) ;
		Fx_ct[i][j][BY] = 0.125*( 
			2.*Fx[i][j][BY]
			+ Fx[i][j+1][BY]
			+ Fx[i][j-1][BY]
			- Fy[i][j][BX]
			- Fy[i-1][j][BX]
			- Fy[i][j+1][BX]
			- Fy[i-1][j+1][BX] ) ;
		Fy_ct[i][j][BY] = 0. ;
	}
	LOOP {
		Fx[i][j][BX] = Fx_ct[i][j][BX] ;
		Fy[i][j][BX] = Fy_ct[i][j][BX] ;
		Fx[i][j][BY] = Fx_ct[i][j][BY] ;
		Fy[i][j][BY] = Fy_ct[i][j][BY] ;
	}
	bounds() ;
#endif

	LOOP {
		/* calculate conserved quantities */
		primtoU(p[i][j], U) ;

		/* evolve conserved quantities */
		PLOOP   U[k] -= (dt/dx)*(Fx[i+1][j][k] - Fx[i][j][k]) 
			+ (dt/dy)*(Fy[i][j+1][k] - Fy[i][j][k]) ;

		/* recover new primitive variables */
		Utoprim(U, p[i][j]) ;

	}
	/* apply boundary conditions */
	bounds() ;

	/* done! */

	t += dt ;
}

