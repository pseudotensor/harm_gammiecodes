#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void usrfun();
#define FREERETURNBAD {free_dmatrix(fjac,1,n,1,n);free_dvector(fvec,1,n);\
	free_dvector(p,1,n);free_ivector(indx,1,n);return(0);}
#define FREERETURNGOOD {free_dmatrix(fjac,1,n,1,n);free_dvector(fvec,1,n);\
	free_dvector(p,1,n);free_ivector(indx,1,n);return(1);}

int mnewt(ntrial,x,n,tolx,tolf)
double tolf,tolx,x[];
int n,ntrial;
{
	void lubksb(),ludcmp();
	int k,i,*indx;
	double errx,errf,d,*fvec,**fjac,*p;

	indx=ivector(1,n);
	p=dvector(1,n);
	fvec=dvector(1,n);
	fjac=dmatrix(1,n,1,n);
	for (k=1;k<=ntrial;k++) {
		/*
		fprintf(stderr,"%d ",k) ;
		*/
		usrfun(x,n,fvec,fjac);
		errf=0.0;
		for (i=1;i<=n;i++) errf += fabs(fvec[i]);
		if (errf <= tolf) {
			/*
			fprintf(stderr,"%d \n",k) ;
			*/
			FREERETURNGOOD
		}
		for (i=1;i<=n;i++) p[i] = -fvec[i];
		ludcmp(fjac,n,indx,&d);
		lubksb(fjac,n,indx,p);
		errx=0.0;
		for (i=1;i<=n;i++) {
			errx += fabs(p[i]);
			x[i] += p[i];
		}
		if (errx <= tolx) {
			/*
			fprintf(stderr,"%d \n",k) ;
			*/
			FREERETURNGOOD
		}
	}
	fprintf(stderr,"mnewt: %d %g\n",k,x[1]) ;
	FREERETURNBAD
}
#undef FREERETURN
