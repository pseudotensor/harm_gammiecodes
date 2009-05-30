
/* 
   Do transport step.  Order of sweep direction is
   varied from step to step. 
*/

#include "decs.h"

void step_trans()
{
	static int nstep = 0 ;
	void sweepx() ;
	void sweepy() ;

	if(nstep%2 == 0) {
		sweepx() ;
		sweepy() ;
	}
	else {
		sweepy() ;
		sweepx() ;
	}

	nstep++ ;
}
