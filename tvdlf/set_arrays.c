
#include "decs.h"

void set_arrays()
{

	p =   (double (*) [NY+2][NP])(& (  a_p[1][1][0])) ;
	ph =  (double (*) [NY+2][NP])(& ( a_ph[1][1][0])) ;
	Fx =  (double (*) [NY+2][NP])(& ( a_Fx[1][1][0])) ;
	Fy =  (double (*) [NY+2][NP])(& ( a_Fy[1][1][0])) ;
	dqx = (double (*) [NY+2][NP])(& (a_dqx[1][1][0])) ;
	dqy = (double (*) [NY+2][NP])(& (a_dqy[1][1][0])) ;

}
