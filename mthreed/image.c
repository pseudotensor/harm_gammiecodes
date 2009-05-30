/*
	produces a ppm file.  
*/

#include "decs.h"

void image(FILE *fp)
{
	double iq,liq,a,lmax,lmin,avg,rms ;
	void john_pal(double data, double min, double max,
                   int *pRed, int *pGreen, int *pBlue);
	int i,j,k ;
	int red,green,blue ;

	/* slice choice */
	k = NZ/2 ;

	/* ppm header */
	fprintf(fp,"P6\n%d %d\n255\n",NX,NY);

	/* density mapping is linear between max dens and min */
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) {
		iq = rho[i][j][k] ;
		avg += iq ;
		rms += iq*iq ;
	}
	avg /= NX*NY ;
	rms /= NX*NY ;
	rms += 1.e-4 ;
	rms = sqrt(rms - avg*avg) ;

	fprintf(stderr,"avg,rms: %g %g\n",avg,rms) ;

	lmax = avg + 3.*rms ;
	lmin = avg - 3.*rms ;

	a = 255./(lmax - lmin) ;
	for(j=NY-1;j>=0;j--) 
	for(i=0;i<NX;i++) {
		iq = rho[i][j][k] ;
		liq = a*(iq - lmin) ;
		if(liq > 255.) liq = 255. ;
		if(liq < 0.) liq = 0. ;

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
