
#include "decs.h"

void stepvar()
{
	void aupdate() ;
	void moc_ct() ;
	void potsolv() ;
	void step_bz() ;
	void step_pg() ;
	void step_visc() ;
	void step_grid() ;
	void step_ie() ;
	void step_trans() ;
	void step_res() ;

	/* grid variation */
	step_grid() ;

	/* solve for potential */
	if(G > GMIN) potsolv(rho,pot) ;

	/* source steps */
	step_pg() ; 
	step_visc() ;
	if(fabs(gam - 1.) > 1.e-6) step_ie() ; 
	step_bz() ; 
	step_res() ;

	moc_ct() ;

	/* forcing */
	if(t+dt >= tdrive && ampf > 0.) {
		step_drive() ;
		tdrive += dtdrive ;
	}

	/* transport steps */
	step_trans() ; 

	t += dt ;
	aupdate() ;
}
