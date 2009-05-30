
/* step magnetic field, velocity using MOC-CT */

#include "decs.h"

void moc_ct()
{
	void bsweepx(),bsweepy(),bsweepz() ;

	/* do field updates by sweeping through planes.
	   interchange order of sweeps */
	/*
	bsweepx() ; 
	bsweepy() ; 
	bsweepz() ; 
	*/
	if(nstep%6 == 0) { bsweepx() ; bsweepy() ; bsweepz() ; }
	if(nstep%6 == 1) { bsweepy() ; bsweepz() ; bsweepx() ; }
	if(nstep%6 == 2) { bsweepz() ; bsweepx() ; bsweepy() ; }
	if(nstep%6 == 3) { bsweepy() ; bsweepx() ; bsweepz() ; }
	if(nstep%6 == 4) { bsweepx() ; bsweepz() ; bsweepy() ; }
	if(nstep%6 == 5) { bsweepz() ; bsweepy() ; bsweepx() ; }
}

void bsweepz()
{
	void bound_var( REAL (*var)[NY+5][NZ+5], REAL t, int sym) ;
	void emf_calc( REAL (*v1)[NN+4], REAL (*v2)[NN+4], REAL (*b1)[NN+4], REAL (*b2)[NN+4], 
		REAL (*r)[NN+4], REAL (*emf2)[NN+1], REAL (*b1prim2)[NN+1], REAL (*b2prim2)[NN+1], 
		int n1, int n2, REAL d1, REAL d2)  ;
	static REAL v1as[NN+4][NN+4], v2as[NN+4][NN+4], 
		b1as[NN+4][NN+4], b2as[NN+4][NN+4], ras[NN+4][NN+4],
		b1prim2as[NN+1][NN+1], b2prim2as[NN+1][NN+1], emf2as[NN+1][NN+1] ;
	static REAL (*v1)[NN+4], (*v2)[NN+4], (*b1)[NN+4], (*b2)[NN+4], (*r)[NN+4] ;
	static REAL (*b1prim2)[NN+1], (*b2prim2)[NN+1], (*emf2)[NN+1] ;
	REAL (*emf)[NY+4][NZ+4],(*b1prim)[NY+4][NZ+4],(*b2prim)[NY+4][NZ+4] ;
	REAL ra,bxa,bya,bza ;
	int i,j,k ;

	/* pointer offset */
	v1 =      (REAL (*)[NN+4])( &(v1as[2][2])) ;
	v2 =      (REAL (*)[NN+4])( &(v2as[2][2])) ;
	b1 =      (REAL (*)[NN+4])( &(b1as[2][2])) ;
	b2 =      (REAL (*)[NN+4])( &(b2as[2][2])) ;
	b1prim2 = (REAL (*)[NN+1])( &(b1prim2as[0][0])) ;
	b2prim2 = (REAL (*)[NN+1])( &(b2prim2as[0][0])) ;
	emf2 =    (REAL (*)[NN+1])( &(emf2as[0][0])) ;
	r =       (REAL (*)[NN+4])( &(ras[2][2])) ;

	emf = work1 ;
	b1prim = work2 ;
	b2prim = work3 ;

	/* x,y plane */
	for(k=0;k<NZ;k++) {
		/* copy plane into 2D array */
		for(i=-2;i<NX+2;i++)
		for(j=-2;j<NY+2;j++) {
			v1[i][j] = vx[i][j][k] ;
			v2[i][j] = vy[i][j][k] ;
			b1[i][j] = bx[i][j][k] ;
			b2[i][j] = by[i][j][k] ;
			r[i][j] = rho[i][j][k] ;
		}
		emf_calc(v1,v2,b1,b2,r,emf2,b1prim2,b2prim2,NX,NY,dx,dy) ;
		for(i=0;i<NX+1;i++)
		for(j=0;j<NY+1;j++) {
			emf[i][j][k] = emf2[i][j] ;
			b1prim[i][j][k] = b1prim2[i][j] ;
			b2prim[i][j][k] = b2prim2[i][j] ;
		}
	}

	/* evolve active fields, velocities */
	/* update field */
	LOOP {
		bx[i][j][k] += dt*(emf[i][j+1][k] - emf[i][j][k])/dy  ;
		by[i][j][k] += dt*(-emf[i+1][j][k] + emf[i][j][k])/dx  ;
	}

	/* update velocity */
	LOOP {
		ra = 0.5*(rho[i][j][k]+rho[i-1][j][k]) ;
		bya = 0.25*(by[i][j][k] + by[i-1][j][k] + 
			by[i][j+1][k] + by[i-1][j+1][k]) ;

		vx[i][j][k] += dt*bya*(b1prim[i][j+1][k] - b1prim[i][j][k])/(dy*ra) ;

		ra = 0.5*(rho[i][j][k]+rho[i][j-1][k]) ;
		bxa = 0.25*(bx[i][j][k] + bx[i+1][j][k] + 
			bx[i][j-1][k] + bx[i+1][j-1][k]) ;

		vy[i][j][k] += dt*bxa*(b2prim[i+1][j][k] - b2prim[i][j][k])/(dx*ra) ;
	}

	/* now: remap magnetic & velocity field */
	bound_var(vx,t,SYMMETRIC) ;
	bound_var(vy,t,SYMMETRIC) ;
	bound_var(by,t,SYMMETRIC) ;
        bound_var(bx,t,SYMMETRIC) ;
}

void bsweepx()
{
	void bound_var( REAL (*var)[NY+5][NZ+5], REAL t, int sym) ;
	void emf_calc( REAL (*v1)[NN+4], REAL (*v2)[NN+4], REAL (*b1)[NN+4], REAL (*b2)[NN+4], 
		REAL (*r)[NN+4], REAL (*emf2)[NN+1], REAL (*b1prim2)[NN+1], REAL (*b2prim2)[NN+1], 
		int n1, int n2, REAL d1, REAL d2)  ;
	static REAL v1as[NN+4][NN+4], v2as[NN+4][NN+4], 
		b1as[NN+4][NN+4], b2as[NN+4][NN+4], ras[NN+4][NN+4],
		b1prim2as[NN+1][NN+1], b2prim2as[NN+1][NN+1], emf2as[NN+1][NN+1] ;
	static REAL (*v1)[NN+4], (*v2)[NN+4], (*b1)[NN+4], (*b2)[NN+4], (*r)[NN+4] ;
	static REAL (*b1prim2)[NN+1], (*b2prim2)[NN+1], (*emf2)[NN+1] ;
	REAL (*emf)[NY+4][NZ+4],(*b1prim)[NY+4][NZ+4],(*b2prim)[NY+4][NZ+4] ;
	REAL ra,bxa,bya,bza ;
	int i,j,k ;

	/* pointer offset */
	v1 =      (REAL (*)[NN+4])( &(v1as[2][2])) ;
	v2 =      (REAL (*)[NN+4])( &(v2as[2][2])) ;
	b1 =      (REAL (*)[NN+4])( &(b1as[2][2])) ;
	b2 =      (REAL (*)[NN+4])( &(b2as[2][2])) ;
	b1prim2 = (REAL (*)[NN+1])( &(b1prim2as[0][0])) ;
	b2prim2 = (REAL (*)[NN+1])( &(b2prim2as[0][0])) ;
	emf2 =    (REAL (*)[NN+1])( &(emf2as[0][0])) ;
	r =       (REAL (*)[NN+4])( &(ras[2][2])) ;

	emf = work1 ;
	b1prim = work2 ;
	b2prim = work3 ;

	/* y,z plane */
	for(i=0;i<NX;i++) {
		/* copy plane into 2D array */
		for(j=-2;j<NY+2;j++) 
		for(k=-2;k<NZ+2;k++) {
			v1[j][k] = vy[i][j][k] ;
			v2[j][k] = vz[i][j][k] ;
			b1[j][k] = by[i][j][k] ;
			b2[j][k] = bz[i][j][k] ;
			r[j][k] = rho[i][j][k] ;
		}
		emf_calc(v1,v2,b1,b2,r,emf2,b1prim2,b2prim2,NY,NZ,dy,dz) ;
		for(j=0;j<NY+1;j++) 
		for(k=0;k<NZ+1;k++) {
			emf[i][j][k] = emf2[j][k] ;
			b1prim[i][j][k] = b1prim2[j][k] ;
			b2prim[i][j][k] = b2prim2[j][k] ;
		}
	}

	/* evolve active fields, velocities */
	/* update field */
	LOOP {
		by[i][j][k] += dt*(emf[i][j][k+1] - emf[i][j][k])/dz  ;
		bz[i][j][k] += dt*(-emf[i][j+1][k] + emf[i][j][k])/dy  ;
	}

	/* update velocity */
	LOOP {
		ra = 0.5*(rho[i][j][k]+rho[i][j-1][k]) ;
		bza = 0.25*(bz[i][j][k] + bz[i][j-1][k] + 
			bz[i][j][k+1] + bz[i][j-1][k+1]) ;

		vy[i][j][k] += dt*bza*(b1prim[i][j][k+1] - b1prim[i][j][k])/(dz*ra) ;

		ra = 0.5*(rho[i][j][k]+rho[i][j][k-1]) ;
		bya = 0.25*(by[i][j][k] + by[i][j+1][k] + 
			by[i][j][k-1] + by[i][j+1][k-1]) ;

		vz[i][j][k] += dt*bya*(b2prim[i][j+1][k] - b2prim[i][j][k])/(dy*ra) ;
	}

	/* now: remap magnetic & velocity field */
	bound_var(vy,t,SYMMETRIC) ;
	bound_var(vz,t,ANTISYMMETRIC) ;
	bound_var(by,t,SYMMETRIC) ;
	bound_var(bz,t,ANTISYMMETRIC) ;
}

void bsweepy()
{
	void bound_var( REAL (*var)[NY+5][NZ+5], REAL t, int sym) ;
	void emf_calc( REAL (*v1)[NN+4], REAL (*v2)[NN+4], REAL (*b1)[NN+4], REAL (*b2)[NN+4], 
		REAL (*r)[NN+4], REAL (*emf2)[NN+1], REAL (*b1prim2)[NN+1], REAL (*b2prim2)[NN+1], 
		int n1, int n2, REAL d1, REAL d2)  ;
	static REAL v1as[NN+4][NN+4], v2as[NN+4][NN+4], 
		b1as[NN+4][NN+4], b2as[NN+4][NN+4], ras[NN+4][NN+4],
		b1prim2as[NN+1][NN+1], b2prim2as[NN+1][NN+1], emf2as[NN+1][NN+1] ;
	static REAL (*v1)[NN+4], (*v2)[NN+4], (*b1)[NN+4], (*b2)[NN+4], (*r)[NN+4] ;
	static REAL (*b1prim2)[NN+1], (*b2prim2)[NN+1], (*emf2)[NN+1] ;
	REAL (*emf)[NY+4][NZ+4],(*b1prim)[NY+4][NZ+4],(*b2prim)[NY+4][NZ+4] ;
	REAL ra,bxa,bya,bza ;
	int i,j,k ;

	/* pointer offset */
	v1 =      (REAL (*)[NN+4])( &(v1as[2][2])) ;
	v2 =      (REAL (*)[NN+4])( &(v2as[2][2])) ;
	b1 =      (REAL (*)[NN+4])( &(b1as[2][2])) ;
	b2 =      (REAL (*)[NN+4])( &(b2as[2][2])) ;
	r =       (REAL (*)[NN+4])( &(ras[2][2])) ;
	b1prim2 = (REAL (*)[NN+1])( &(b1prim2as[0][0])) ;
	b2prim2 = (REAL (*)[NN+1])( &(b2prim2as[0][0])) ;
	emf2 =    (REAL (*)[NN+1])( &(emf2as[0][0])) ;

	emf = work1 ;
	b1prim = work2 ;
	b2prim = work3 ;

	/* x,z plane */
	for(j=0;j<NY;j++) {
		/* copy plane into 2D array */
		for(i=-2;i<NX+2;i++) 
		for(k=-2;k<NZ+2;k++) {
			v1[i][k] = vx[i][j][k] ;
			v2[i][k] = vz[i][j][k] ;
			b1[i][k] = bx[i][j][k] ;
			b2[i][k] = bz[i][j][k] ;
			r[i][k] = rho[i][j][k] ;
		}
		emf_calc(v1,v2,b1,b2,r,emf2,b1prim2,b2prim2,NX,NZ,dx,dz) ;
		for(i=0;i<NX+1;i++) 
		for(k=0;k<NZ+1;k++) {
			emf[i][j][k] = emf2[i][k] ;
			b1prim[i][j][k] = b1prim2[i][k] ;
			b2prim[i][j][k] = b2prim2[i][k] ;
		}
	}

	/* evolve active fields, velocities */
	/* update field */
	LOOP {
		bx[i][j][k] += dt*(emf[i][j][k+1] - emf[i][j][k])/dz  ;
		bz[i][j][k] += dt*(-emf[i+1][j][k] + emf[i][j][k])/dx  ;
	}

	/* update velocity */
	LOOP {
		ra = 0.5*(rho[i][j][k]+rho[i-1][j][k]) ;
		bza = 0.25*(bz[i][j][k] + bz[i-1][j][k] + 
			bz[i][j][k+1] + bz[i-1][j][k+1]) ;

		vx[i][j][k] += dt*bza*(b1prim[i][j][k+1] - b1prim[i][j][k])/(dz*ra) ;

		ra = 0.5*(rho[i][j][k]+rho[i][j][k-1]) ;
		bxa = 0.25*(bx[i][j][k] + bx[i+1][j][k] + 
			bx[i][j][k-1] + bx[i+1][j][k-1]) ;

		vz[i][j][k] += dt*bxa*(b2prim[i+1][j][k] - b2prim[i][j][k])/(dx*ra) ;
	}

	/* now: remap magnetic & velocity field */
	bound_var(vx,t,SYMMETRIC) ;
	bound_var(vz,t,ANTISYMMETRIC) ;
	bound_var(bz,t,ANTISYMMETRIC) ;
        bound_var(bx,t,SYMMETRIC) ;
}

void emf_calc( REAL (*v1)[NN+4], REAL (*v2)[NN+4], REAL (*b1)[NN+4], REAL (*b2)[NN+4], 
	REAL (*r)[NN+4], REAL (*emf)[NN+1], REAL (*b1prim)[NN+1], REAL (*b2prim)[NN+1], 
	int n1, int n2, REAL d1, REAL d2) 
{
	static REAL dqvas[NN+2][NN+2], dqbas[NN+2][NN+2] ;
	static REAL (*dqv)[NN+2], (*dqb)[NN+2] ;
	static REAL v1star[NN+1][NN+1], v2star[NN+1][NN+1] ;
	static REAL b1star[NN+1][NN+1], b2star[NN+1][NN+1] ;
	REAL b1a,b2a,Dx1,Dx2,Dy1,Dy2 ;
	REAL ra,vaa,v1a,v2a,sra,vp,vm,bp,bm ;
	REAL sgn_va ;
	int i,j ;
	void dqx_calc2(REAL (*var)[NN+4], REAL (*dq)[NN+2], int n1, int n2) ;
	void dqy_calc2(REAL (*var)[NN+4], REAL (*dq)[NN+2], int n1, int n2) ;

	/* pointer shift */
	dqv = (REAL (*)[NN+2])(& (dqvas[1][1])) ;
	dqb = (REAL (*)[NN+2])(& (dqbas[1][1])) ;

	/* sweep in 1-direction */
	/* get slopes */
	dqx_calc2(v2,dqv,n1,n2) ;
	dqx_calc2(b2,dqb,n1,n2) ;

	for(i=0;i<n1+1;i++)
	for(j=0;j<n2+1;j++) {
		ra = 0.25*(r[i][j] + r[i-1][j] +
			r[i][j-1] + r[i-1][j-1]) ;
		sra = sqrt(ra) ;
		vaa = 0.5*((b1[i][j]+b1[i][j-1])/sra) ;
		sgn_va = copysign(1.,vaa) ;
		vaa = fabs(vaa) ;
		v1a = 0.5*(v1[i][j]+v1[i][j-1]) ;
		Dx1 = dt*vaa/d1 ;
		Dx2 = dt*v1a/d1 ;

		if(Dx1 + Dx2 > 0.) {
			vp = v2[i-1][j] + (1. - Dx1 - Dx2)*dqv[i-1][j] ;
			bp = b2[i-1][j] + (1. - Dx1 - Dx2)*dqb[i-1][j] ;
		}
		else {
			vp = v2[i][j] + (-1. - Dx1 - Dx2)*dqv[i][j] ;
			bp = b2[i][j] + (-1. - Dx1 - Dx2)*dqb[i][j] ;
		}
		if(-Dx1 + Dx2 > 0.) {
			vm = v2[i-1][j] + (1. + Dx1 - Dx2)*dqv[i-1][j] ;
			bm = b2[i-1][j] + (1. + Dx1 - Dx2)*dqb[i-1][j] ;
		}
		else {
			vm = v2[i][j] + (-1. + Dx1 - Dx2)*dqv[i][j] ;
			bm = b2[i][j] + (-1. + Dx1 - Dx2)*dqb[i][j] ;
		}
		
                /* solution to Stone & Norman eqtn 43,44-- use constant r */
                v2star[i][j] = 0.5*(vm + vp + sgn_va*(bm - bp)/sra) ;
                b2star[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*sra) ;

		/* find bprim, for updating velocities */
                vp = v2[i-1][j] + (1. - Dx1)*dqv[i-1][j] ;
                bp = b2[i-1][j] + (1. - Dx1)*dqb[i-1][j] ;
                vm = v2[i][j] + (-1. + Dx1)*dqv[i][j] ;
                bm = b2[i][j] + (-1. + Dx1)*dqb[i][j] ;

                b2prim[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*sra) ;
	}

	/* sweep in 2-direction */
	/* first get slopes */
	dqy_calc2(v1,dqv,n1,n2) ;
	dqy_calc2(b1,dqb,n1,n2) ;

	for(i=0;i<n1+1;i++)
	for(j=0;j<n2+1;j++) {
		ra = 0.25*(r[i][j] + r[i-1][j] +
			r[i][j-1] + r[i-1][j-1]) ;
		sra = sqrt(ra) ;
		vaa = 0.5*((b2[i][j]+b2[i-1][j])/sra) ;
		sgn_va = copysign(1.,vaa) ;
		vaa = fabs(vaa) ;
		v2a = 0.5*(v2[i][j]+v2[i-1][j]) ;
		Dy1 = dt*vaa/d2 ;
		Dy2 = dt*v2a/d2 ;

                /* values at the foot of the plus characteristic
                  are, by convention, in i-1 zone */
		if(Dy1 + Dy2 > 0.) {
			vp = v1[i][j-1] + (1. - Dy1 - Dy2)*dqv[i][j-1] ;
			bp = b1[i][j-1] + (1. - Dy1 - Dy2)*dqb[i][j-1] ;
		}
		else {
			vp = v1[i][j] + (-1. - Dy1 - Dy2)*dqv[i][j] ;
			bp = b1[i][j] + (-1. - Dy1 - Dy2)*dqb[i][j] ;
		}
		if(-Dy1 + Dy2 > 0.) {
			vm = v1[i][j-1] + (1. + Dy1 - Dy2)*dqv[i][j-1] ;
			bm = b1[i][j-1] + (1. + Dy1 - Dy2)*dqb[i][j-1] ;
		}
		else {
			vm = v1[i][j] + (-1. + Dy1 - Dy2)*dqv[i][j] ;
			bm = b1[i][j] + (-1. + Dy1 - Dy2)*dqb[i][j] ;
		}
		
                /* solution to Stone & Norman eqtn 43,44-- use constant r */
                v1star[i][j] = 0.5*(vm + vp + sgn_va*(bm - bp)/sra) ;
                b1star[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*sra) ;

		/* find bprim, to update velocities */
                vp = v1[i][j-1] + (1. - Dy1)*dqv[i][j-1] ;
                vm = v1[i][j] + (-1. + Dy1)*dqv[i][j] ;
                bp = b1[i][j-1] + (1. - Dy1)*dqb[i][j-1] ;
                bm = b1[i][j] + (-1. + Dy1)*dqb[i][j] ;

                b1prim[i][j] = 0.5*(bm + bp + sgn_va*(vm - vp)*sra) ;
	}

	/* calculate emf */
	for(i=0;i<n1+1;i++)
	for(j=0;j<n2+1;j++) {
		emf[i][j] = v1star[i][j]*b2star[i][j] - v2star[i][j]*b1star[i][j] ;
	}
}

void dqx_calc2(REAL (*var)[NN+4], REAL (*dq)[NN+2], int n1, int n2) 
{
	int i,j ;
	REAL Dqm,Dqp,pr ;

	for(i=-1;i<n1+1;i++) 
	for(j=-1;j<n2+1;j++) {
		Dqm = var[i][j]-var[i-1][j] ;
		Dqp = var[i+1][j]-var[i][j] ;

		pr = Dqm*Dqp ;
		dq[i][j] = (pr > 0.) ? pr/(Dqm+Dqp) : 0. ;
	}
}

void dqy_calc2(REAL (*var)[NN+4], REAL (*dq)[NN+2], int n1, int n2) 
{
	int i,j ;
	REAL Dqm,Dqp,pr ;

	for(i=-1;i<n1+1;i++) 
	for(j=-1;j<n2+1;j++) {
		Dqm = var[i][j]-var[i][j-1] ;
		Dqp = var[i][j+1]-var[i][j] ;

		pr = Dqm*Dqp ;
		dq[i][j] = (pr > 0.) ? pr/(Dqm+Dqp) : 0. ;
	}
}
