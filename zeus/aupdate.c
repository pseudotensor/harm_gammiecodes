
/* update expansion factor a */

#include "decs.h"

void aupdate()
{
	
	if(tdecay < 1.e6) {
		a = exp(-t/tdecay) ;
		da = -a/tdecay ;
	}
	else {
		a = 1. ;
		da = 0. ;
	}
}
