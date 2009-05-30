
#include "decs.h"
#include "defs.h"

main(int argc,char *argv[])
{
	void init(),timestep(),stepvar(),dump(),diag(),potsolv() ;
	double tnext ;
	static int nstep = 0 ;

	tnext = 0. ;

	sscanf(argv[1],"%d",&NX) ;
	NY = NX ;

	fprintf(stderr,"starting...\n") ;

	/* perform initializations */
	init() ;

	fprintf(stderr,"done init..\n") ;

	/* find potential */
	if(G > GMIN) potsolv(rho,pot) ;
	
	fprintf(stderr,"done pot..\n") ;

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

		fprintf(stderr,"t,dt: %10.5g %10.5g\n",t,dt) ;
		nstep++ ;
	}

	/* dump state to stdout */
	/*
	dump(stdout) ;
	*/

	/* do final diagnostics */
	diag(2) ;

	fprintf(stderr,"nstep, nz: %d %d\n",nstep,nstep*NX*NY) ;

	return(0) ;
}

