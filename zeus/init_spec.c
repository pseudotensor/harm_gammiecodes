
#include "decs.h"
#define KMAX	8

void init_spec()
{
	void gauss_gen() ;
	double norm,**dmatrix();
	double **ax,**ay,**az ;
	double n ;
	int i,j,ip,jp ;

	ax = dmatrix(0,NX-1,0,NY-1) ;
	ay = dmatrix(0,NX-1,0,NY-1) ;
	az = dmatrix(0,NX-1,0,NY-1) ;

	/* generate random fields */
	n = 2. ;	/* E_k \sim k^(-n) */
	gauss_gen(ax,-0.5*(n + 3.)) ;
	gauss_gen(ay,-0.5*(n + 3.)) ;
	gauss_gen(az,-0.5*(n + 3.)) ;

	/* difference potential to get divergence-free velocity */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		ip = i+1 ;
		jp = j+1 ;
		if(ip == NX) ip = 0 ;
		if(jp == NY) jp = 0 ;
		vx[i][j] = (az[i][jp] - az[i][j])/(a*dy) ;
		vy[i][j] = -(az[ip][j] - az[i][j])/(a*dx) ;
		vz[i][j] = (ay[ip][j] - ay[i][j])/(a*dx) +
			(ax[i][jp] - ax[i][j])/(a*dy) ;
	}

	/* do same for the magnetic field */
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		Bx[i][j] = 0. ;
		By[i][j] = 0. ;
		Bz[i][j] = 0. ;
	}
#if 0
	gauss_gen(ax,-2.5) ;
	gauss_gen(ay,-2.5) ;
	gauss_gen(az,-2.5) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		ip = i+1 ;
		jp = j+1 ;
		if(ip == NX) ip = 0 ;
		if(jp == NY) jp = 0 ;
		Bx[i][j] = (az[i][jp] - az[i][j])/(a*dy) ;
		By[i][j] = -(az[ip][j] - az[i][j])/(a*dx) ;
		Bz[i][j] = (ay[ip][j] - ay[i][j])/(a*dx) +
			(ax[i][jp] - ax[i][j])/(a*dy) ;
	}
#endif

	/* normalize */
	norm = 0. ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		/* assumes density = 1 */
		norm += 0.5*(
			vx[i][j]*vx[i][j] +
			vy[i][j]*vy[i][j] +
			vz[i][j]*vz[i][j]) +
			0.5*(
			Bx[i][j]*Bx[i][j] +
			By[i][j]*By[i][j] +
			Bz[i][j]*Bz[i][j]) ;
	}
	norm = sqrt(fabs(ampf)/(norm*dx*dy*a*a)) ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vx[i][j] *= norm ;
		vy[i][j] *= norm ;
		vz[i][j] *= norm ;
		Bx[i][j] *= norm ;
		By[i][j] *= norm ;
		Bz[i][j] *= norm ;
	}

}

/* add Gaussian random field to matrix f, with power spectrum
   amplitude norm */
void gauss_gen(f,sl)
double **f ;
double sl ; 
{
	float **dfr,**dfi ;
	float **matrix() ;
	double kx,ky,phase,amp,powsp(),ranc() ;
	void fft_grid() ;
	int i,j,ip,jp ;

	dfr = matrix(0,NX-1,0,NY-1) ;
	dfi = matrix(0,NX-1,0,NY-1) ;

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
			amp = powsp((double)ip,(double)jp,sl) ;
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

double powsp(kx,ky,sl)
double kx,ky,sl ;
{
	double k2,p ;
	double ranc(),x1,x2 ;


	k2 = kx*kx + ky*ky ;
	p = pow(k2,0.5*sl) ;

	x1 = ranc(0) ;
	x2 = ranc(0) ;	/* wasteful */
	p *= sqrt(-2.*log(x1))*cos(2.*M_PI*x2) ;


	return(p) ;
}
