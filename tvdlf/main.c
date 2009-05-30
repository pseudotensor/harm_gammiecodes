
#include "decs.h"
#include "defs.h"

main(argc,argv)
int argc ;
char *argv[] ;
{
	void init(),timestep(),step_ch(),dump(),diag() ;
	double tnext ;
	int nstep ;

	/* warn user about boundary conditions */
	if(PER) fprintf(stderr,"Using PER boundary conditions\n") ;
	else if(STD) fprintf(stderr,"Using STD boundary conditions\n") ;
	else {
		fprintf(stderr,"Unknown boundary conditions\n") ;
		return(1) ;
	}

	/* perform initializations */
	init() ;

	/* do initial diagnostics */
	diag(0) ;

	nstep = 0 ;
	tnext = t ;
	while(t < tf) {

		fprintf(stderr,"%d %g %g\n",nstep,t,dt) ;
		fflush(stderr);
		/* step variables forward in time */
		step_ch() ;
		
		/* perform diagnostics */
		if(t >= tnext) {
			diag(1) ;
			tnext += DTl ;
		}
		nstep++ ;
	}
	fprintf(stderr,"ns,ts: %d %d\n",nstep,nstep*NX*NY) ;

	/* do final diagnostics */
	diag(2) ;
	exit(0) ;
}

