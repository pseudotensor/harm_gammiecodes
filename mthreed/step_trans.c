
/* 
   Do transport step.  Order of sweep direction is
   varied from step to step. 
*/

#include "decs.h"

void step_trans()
{
	void sweepx() ;
	void sweepy() ;
	void sweepym() ;
	void sweepz() ;

	if(nstep%6 == 0) { sweepx() ; sweepy() ; sweepz() ; }
	if(nstep%6 == 1) { sweepz() ; sweepx() ; sweepy() ; }
	if(nstep%6 == 2) { sweepy() ; sweepz() ; sweepx() ; }
	if(nstep%6 == 3) { sweepx() ; sweepz() ; sweepy() ; }
	if(nstep%6 == 4) { sweepy() ; sweepx() ; sweepz() ; }
	if(nstep%6 == 5) { sweepz() ; sweepy() ; sweepx() ; }
#if SWEEPYM
	if(NY > 1 && q!=0.) sweepym() ;
#endif
}
