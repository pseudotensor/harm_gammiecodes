/*
	produces a ppm file.  
*/

#include "decs.h"

void image(FILE *fp)
{
	float iq,liq,a,b,lmax,lmin,X ;
	void john_pal(double data, double min, double max,
                   int *pRed, int *pGreen, int *pBlue);
	int i,k,nx1,nx2 ;
	int red,green,blue ;
	static int slice = NY/2 ;
	static int nghost = 0 ;
	static int nx = NX ;
	static int nz = NZ ;

	nx1 = nx + nghost ;
	nx2 = nz + nghost ;
	
	fprintf(fp,"P6\n%d %d\n255\n",nx1,nx2);

	/* density mapping is logarithmic, in 255 steps
	   between e^lmax and e^lmin */
#if LINEAR || QEFF
	lmax = pert ;
	lmin = -pert ;
#elif RHOCOL || STRIPE
	lmax = 3. ;
	lmin = 1. ;
#elif FLOOP
	lmax = 0.08 ;
	lmin = -0.04 ;
#elif BY0
	lmax = 0. ;
	lmin = -q*W*B0*tf ;
#endif
	for(k=nx2-1;k>=-nghost;k--) 
	for(i=-nghost;i<=nx1-1;i++) {
		X = (i + .5)*dx - .5 ;
#if LINEAR
		iq = (bx[i+1][slice][k]-bx[i][slice][k])/dx
		   + (by[i][slice+1][k]-by[i][slice][k])/dy
		   + (bz[i][slice][k+1]-bz[i][slice][k])/dz ;
		if(iq!=0.) iq = log(fabs(iq)) ;
		else iq = lmin ;
		//iq = vz[i][slice][k] ;
		//iq = vx[i][slicej][k] ;
		//iq = vy[i][slice][k] ;
		//iq = rho[i][slice][k] - 1. ;
#elif RHOCOL || STRIPE
		iq = rho[i][slice][k] ;
		//iq = vx[i][slice][k] ;
#elif FLOOP
		iq = (by[i][slice][k] - by[i-1][slice][k])/dx - (bx[i][slice][k] - bx[i][slice-1][k])/dy ;
#elif BY0
		iq = by[i][slice][k] ;
#endif
		john_pal(iq,lmin,lmax,&red,&green,&blue);
      		fputc(red,fp);
      		fputc(green,fp);
      		fputc(blue,fp);
	}
}


void john_pal(double data, double min, double max,
                   int *pRed, int *pGreen, int *pBlue){
  double max_min = max - min;
  double a, b, c, d, e, f;
  double x, y;

  if(max_min > 0.0){ 
    x = (data - min)/(max_min); 
    
    /* ========== Red ============ */
    a = 4.0*x - 1.52549019607844;
    b = 4.52941176470589 - 4.0*x;
    y = a < b ? a : b;
    *pRed = (int)(255.0*y);
    *pRed = *pRed >   0 ? *pRed :   0;
    *pRed = *pRed < 255 ? *pRed : 255;

    /* ========== Green ========== */
    a = 4.0*x - 0.521568627450979;
    b = 2.52549019607844 - 4.0*x;
    c = a < b ? a : b;
    d = 4.0*x - 1.53725490196073;
    e = 3.52941176470581 - 4.0*x;
    f = d < e ? d : e;
    y = c > f ? c : f;
    *pGreen = (int)(255.0*y);
    *pGreen = *pGreen >   0 ? *pGreen :   0;
    *pGreen = *pGreen < 255 ? *pGreen : 255;

    /* ========== Blue =========== */
    a = 4.0*x + 0.498039215686276;
    b = 2.50980392156862 - 4.0*x;
    y = a < b ? a : b;
    *pBlue = (int)(255.0*y);
    *pBlue = *pBlue >   0 ? *pBlue :   0;
    *pBlue = *pBlue < 255 ? *pBlue : 255;
  }
  else{
    *pRed = *pGreen = *pBlue = (data > max ? 255 : 0);
  }

  return;
}
