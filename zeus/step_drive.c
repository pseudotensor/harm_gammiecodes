
#include "decs.h"

/* driving power spectrum cut off at kmax */
#define KMAX	16

void step_drive()
{
	void dgauss_gen() ;
	double norm,**dmatrix();
	static double **ax,**ay,**az ;
	static double **dvx,**dvy,**dvz ;
	double c1,c2,c3x,c3y,c3z,c32,c4 ;
	double alpha,betax,betay,betaz ;
	int i,j,ip,jp ;
	static int firstc = 1 ;

	if(firstc) {
		firstc = 0 ;
		ax = dmatrix(0,NX-1,0,NY-1) ;
		ay = dmatrix(0,NX-1,0,NY-1) ;
		az = dmatrix(0,NX-1,0,NY-1) ;
		dvx = dmatrix(0,NX-1,0,NY-1) ;
		dvy = dmatrix(0,NX-1,0,NY-1) ;
		dvz = dmatrix(0,NX-1,0,NY-1) ;
	}

	/* generate random fields */
	dgauss_gen(ax) ;
	dgauss_gen(ay) ;
	dgauss_gen(az) ;

	/* difference potential to get divergence-free velocity */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		ip = i+1 ;
		jp = j+1 ;
		if(ip == NX) ip = 0 ;
		if(jp == NY) jp = 0 ;
		dvx[i][j] = (az[i][jp] - az[i][j])/(a*dy) ;
		dvy[i][j] = -(az[ip][j] - az[i][j])/(a*dx) ;
		dvz[i][j] = (ay[ip][j] - ay[i][j])/(a*dx) +
			(ax[i][jp] - ax[i][j])/(a*dy) ;
	}

	/* normalize; calculate quadratic coefficients */
	c1 = c2 = c3x = c3y = c3z = c4 = 0. ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		c1 += rho[i][j]*(
			vx[i][j]*dvx[i][j] +
			vy[i][j]*dvy[i][j] +
			vz[i][j]*dvz[i][j]
			) ;
		c2 += 0.5*rho[i][j]*(
			dvx[i][j]*dvx[i][j] +
			dvy[i][j]*dvy[i][j] +
			dvz[i][j]*dvz[i][j]
			) ;
		c3x += rho[i][j]*dvx[i][j] ;
		c3y += rho[i][j]*dvy[i][j] ;
		c3z += rho[i][j]*dvz[i][j] ;

		c4 += rho[i][j] ;
	}
	c1 *= dx*dy ;
	c2 *= dx*dy ;
	c3x *= dx*dy ;
	c3y *= dx*dy ;
	c3z *= dx*dy ;
	c4 *= dx*dy ;
	c32 = c3x*c3x + c3y*c3y + c3z*c3z ;

	/*
	alpha = 0.5*(-c1 + sqrt(c1*c1 + 4.*(c2 - 0.5*c32/c4)*ampf*dtdrive))/
			(c2 - 0.5*c32/c4) ;
	*/
	alpha = (c4*c1 + sqrt(c4*(-2.*ampf*dtdrive*(c32 - 2.*c2*c4)
		+ c4*c1*c1)))/(c32 - 2.*c2*c4) ;
	betax = -alpha*c3x/c4 ;
	betay = -alpha*c3y/c4 ;
	betaz = -alpha*c3z/c4 ;

	/* transform */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vx[i][j] += alpha*dvx[i][j] + betax ;
		vy[i][j] += alpha*dvy[i][j] + betay ;
		vz[i][j] += alpha*dvz[i][j] + betaz ;
	}
}

/* add Gaussian random field to matrix f, with power spectrum
   amplitude norm */
void dgauss_gen(f)
double **f ;
{
	static float **dfr,**dfi ;
	float **matrix() ;
	double kx,ky,phase,amp,dpowsp(),ranc() ;
	void fft_grid() ;
	int i,j,ip,jp ;
	static int firstc = 1 ;

	if(firstc) {
		firstc = 0 ;
		dfr = matrix(0,NX-1,0,NY-1) ;
		dfi = matrix(0,NX-1,0,NY-1) ;
	}

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) dfr[i][j] = dfi[i][j] = 0. ;

	/* generate complex random field; use only real part */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		if(i > NX/2) ip = NX-i ;
		else ip = i ;
		if(j > NY/2) jp = NY-j ;
		else jp = j ;

		if((ip != 0 || jp != 0) && ip*ip+jp*jp < KMAX*KMAX) {
			phase = 2.*M_PI*ranc(0) ;
			amp = dpowsp((double)ip,(double)jp) ;
		}
		else {
			phase = 0. ;
			amp = 0. ;
		}

		dfr[i][j] = (float)amp*cos(phase) ;
		dfi[i][j] = (float)amp*sin(phase) ;
	}

	/* fft */
	fft_grid(dfr,dfi,dfr,dfi,NX,NY,-1) ;

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) f[i][j] = (double)dfr[i][j] ;

	return ;
}

/* this gives E_k ~ k^8 exp(-8 k/kpk) */

#define POWFUNC (k2*exp(-4.*k/kpk))

double dpowsp(kx,ky)
double kx,ky ;
{
	double k2,k,p ;
	double ranc(),x1,x2 ;

	k2 = kx*kx + ky*ky ;
	k = sqrt(k2) ;
	p = POWFUNC ;

	x1 = ranc(0) ;
	x2 = ranc(0) ;	/* wasteful */
	p *= sqrt(-2.*log(x1))*cos(2.*M_PI*x2) ;

	return(p) ;
}
