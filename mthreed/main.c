
#include "decs.h"
#include "defs.h"

main(int argc, char *argv[])
{
	void init(),timestep(),stepvar(),dump(),diag() ;
	REAL tnext ;
	int seed,i,j,k ;

	tnext = 0. ;
	nstep = 0 ;

	if(argc < 1) {
		fprintf(stderr,"usage: bend3d \n") ;
		exit(0) ;
	}

	/* perform initializations */
	init() ;

	/* do initial diagnostics */
	diag(0) ;

	while(t < tf) {

		/* find timestep */
		timestep() ;

		/* step variables forward in time */
		stepvar() ;

		/* perform diagnostics */
		if(t > tnext) {
			diag(1) ;
			tnext += DTl ;
		
		}
		
		//if(nstep == 64) exit(0) ;

		fprintf(stderr,"t,dt: %10.5g %10.5g\n",t,dt) ;
		nstep++ ;
	}

	/* do final diagnostics */
	diag(2) ;

	fprintf(stderr,"nstep: %d\n",nstep) ;
	fprintf(stderr,"zone cycles: %g\n",(double)(NX*NY*NZ)*
		(double)(nstep)) ;

	return(0) ;
}
