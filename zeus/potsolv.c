
#include "decs.h"

/* option to use full Green's function or 
   grid Green's function */
#define	GRIDGREEN	1
#define	FULLGREEN	0

void potsolv(dens,psi)
double **dens,**psi ;
{
	double kx,ky ;
	int i,j,ip,jp ;
	float **matrix() ;
	void fft_grid() ;

	static float **dr,**di,**pr,**pi ;
	static float **kdr,**kdi,**kpr,**kpi ;
	static float **kern ;
	static int first_call = 1 ;

	/* set up arrays */
	if(first_call == 1) {
		first_call = 0 ;
		dr = matrix(0,NX-1,0,NY-1) ;
		di = matrix(0,NX-1,0,NY-1) ;
		pr = matrix(0,NX-1,0,NY-1) ;
		pi = matrix(0,NX-1,0,NY-1) ;
		kdr = matrix(0,NX-1,0,NY-1) ;
		kdi = matrix(0,NX-1,0,NY-1) ;
		kpr = matrix(0,NX-1,0,NY-1) ;
		kpi = matrix(0,NX-1,0,NY-1) ;
		kern = matrix(0,NX-1,0,NY-1) ;

		/* precompute gravitational kernel */
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) {
			ip = (i + NX/2)%NX - NX/2 ;
			jp = (j + NY/2)%NY - NY/2 ;
			if(ip != 0 || jp != 0) {
				kx = 2.*M_PI*ip/(NX*a*dx) ;
				ky = 2.*M_PI*jp/(NY*a*dy) ;

#if GRIDGREEN
				kern[i][j] = -2.*M_PI*G/
					((1. - cos(kx*a*dx))/(a*dx*a*dx) +
					 (1. - cos(ky*a*dy))/(a*dy*a*dy)) ;
#endif
#if FULLGREEN
				kern[i][j] = -4.*M_PI*G/(kx*kx + ky*ky) ;
#endif
			}
			else {
				kern[i][j] = 0. ;
			}
		}
	}

        for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		dr[i][j] = dens[i][j] ;
		di[i][j] = 0.0 ;
	}

	/* forward fft */
	fft_grid(dr,di,kdr,kdi,NX,NY,1) ;

	/* solve Poisson's equation, using Green's
	   function appropriate for grid */
        for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		kpr[i][j] = kern[i][j]*kdr[i][j] ;
		kpi[i][j] = kern[i][j]*kdi[i][j] ;
	}

	/* inverse fft */
	fft_grid(kpr,kpi,pr,pi,NX,NY,-1) ;

        for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		psi[i][j] = (double)pr[i][j] ;
	}

	/* done! */
}

/* fft on grid */
void fft_grid(inr,ini,outr,outi,nx,ny,dir)
float **inr,**ini,**outr,**outi ;
int nx,ny,dir ;
{
	static float *inr1,*ini1,*outr1,*outi1 ;
	static int firstc = 1 ;
	float *vector() ;
	int i,j,n ;
	int tfft() ;

	/* hold onto allocated arrays */
	if(firstc) {
		firstc = 0 ;
		if (nx > ny) n = nx ;
		else n = ny ;

		inr1 = vector(0,n-1) ;
		ini1 = vector(0,n-1) ;
		outr1 = vector(0,n-1) ;
		outi1 = vector(0,n-1) ;
	}

	/* sweep in x direction */
	for(i=0;i<NX;i++) {
		for(j=0;j<NY;j++) {
			inr1[j] = inr[i][j] ;
			ini1[j] = ini[i][j] ;
		}
		tfft(inr1,ini1,outr1,outi1,NY,dir) ;
		for(j=0;j<NY;j++) {
			outr[i][j] = outr1[j];
			outi[i][j] = outi1[j];
		}
	}

	/* sweep in y direction */
	for(j=0;j<NY;j++) {
		for(i=0;i<NX;i++) {
			inr1[i] = outr[i][j] ;
			ini1[i] = outi[i][j] ;
		}
		tfft(inr1,ini1,outr1,outi1,NX,dir) ;
		for(i=0;i<NX;i++) {
			outr[i][j] = outr1[i];
			outi[i][j] = outi1[i];
		}
	}

	/* done! */
}
