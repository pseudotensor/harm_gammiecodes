
#define NFFTMAX	1024

#include <stdio.h>

/* interface to real fft routines */

void tfft(float *in_r,float *in_i,float *out_r,float *out_i,int n,int dir)
/* +1: FFT, -1: IFFT */
{
	double rdum[NFFTMAX] ;
	double idum[NFFTMAX] ;
	int i ;

	for(i=0;i<n;i++) {
		rdum[i] = in_r[i] ;
		idum[i] = in_i[i] ;
	}

	if(dir > 0) fft(n,rdum,idum) ;
	if(dir < 0) ifft(n,rdum,idum) ;

	for(i=0;i<n;i++) {
		out_r[i] = rdum[i] ;
		out_i[i] = idum[i] ;
	}
	if(dir < 0) {
		for(i=0;i<n;i++) {
			out_r[i] /= n ;
			out_i[i] /= n ;
		}
	}

}

