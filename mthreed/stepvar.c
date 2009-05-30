
#include "decs.h"

void stepvar()
{
	void moc_ct() ;
	void step_bz() ;
	void step_pg() ;
	void step_visc() ;
	void step_trans() ;
	void step_ie() ;
	void step_res() ;

	/* source steps */
	step_pg() ; 
	step_visc() ;
	if(fabs(gam - 1.) > 1.e-6) step_ie() ;

	moc_ct() ;

	/* forcing */
	/*
	if (drive > 0.) step_drive() ;
	*/

	/* transport steps */
	step_trans() ; 

	t += dt ;
}
