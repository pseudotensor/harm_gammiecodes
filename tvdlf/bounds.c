
#include "decs.h"

#if PER
void bounds()
{
	int i,j,k ;

	for(i=0;i<NX;i++) {
		PLOOP {
			/* lower boundary */
			p  [i][-1][k] = p  [i][NY-1][k] ;
			ph [i][-1][k] = ph [i][NY-1][k] ;
			Fx [i][-1][k] = Fx [i][NY-1][k] ;
			Fy [i][-1][k] = Fy [i][NY-1][k] ;
			dqx[i][-1][k] = dqx[i][NY-1][k] ;
			dqy[i][-1][k] = dqy[i][NY-1][k] ;

			/* upper boundary */
			p  [i][NY][k] = p  [i][0][k] ;
			ph [i][NY][k] = ph [i][0][k] ;
			Fx [i][NY][k] = Fx [i][0][k] ;
			Fy [i][NY][k] = Fy [i][0][k] ;
			dqx[i][NY][k] = dqx[i][0][k] ;
			dqy[i][NY][k] = dqy[i][0][k] ;
		}
	}
	for(j=0;j<NY;j++) {
		PLOOP {
			/* lower boundary */
			p  [-1][j][k] = p  [NX-1][j][k] ;
			ph [-1][j][k] = ph [NX-1][j][k] ;
			Fx [-1][j][k] = Fx [NX-1][j][k] ;
			Fy [-1][j][k] = Fy [NX-1][j][k] ;
			dqx[-1][j][k] = dqx[NX-1][j][k] ;
			dqy[-1][j][k] = dqy[NX-1][j][k] ;

			/* upper boundary */
			p  [NX][j][k] = p  [0][j][k] ;
			ph [NX][j][k] = ph [0][j][k] ;
			Fx [NX][j][k] = Fx [0][j][k] ;
			Fy [NX][j][k] = Fy [0][j][k] ;
			dqx[NX][j][k] = dqx[0][j][k] ;
			dqy[NX][j][k] = dqy[0][j][k] ;
		}
	}
	/* corners */
	PLOOP {
		p  [-1][-1][k] = p  [NX-1][NY-1][k] ;
		ph [-1][-1][k] = ph [NX-1][NY-1][k] ;
		Fx [-1][-1][k] = Fx [NX-1][NY-1][k] ;
		Fy [-1][-1][k] = Fy [NX-1][NY-1][k] ;
		dqx[-1][-1][k] = dqx[NX-1][NY-1][k] ;
		dqy[-1][-1][k] = dqy[NX-1][NY-1][k] ;

		p  [-1][NY][k] = p  [NX-1][0][k] ;
		ph [-1][NY][k] = ph [NX-1][0][k] ;
		Fx [-1][NY][k] = Fx [NX-1][0][k] ;
		Fy [-1][NY][k] = Fy [NX-1][0][k] ;
		dqx[-1][NY][k] = dqx[NX-1][0][k] ;
		dqy[-1][NY][k] = dqy[NX-1][0][k] ;

		p  [NX][NY][k] = p  [0][0][k] ;
		ph [NX][NY][k] = ph [0][0][k] ;
		Fx [NX][NY][k] = Fx [0][0][k] ;
		Fy [NX][NY][k] = Fy [0][0][k] ;
		dqx[NX][NY][k] = dqx[0][0][k] ;
		dqy[NX][NY][k] = dqy[0][0][k] ;

		p  [NX][-1][k] = p  [0][NY-1][k] ;
		ph [NX][-1][k] = ph [0][NY-1][k] ;
		Fx [NX][-1][k] = Fx [0][NY-1][k] ;
		Fy [NX][-1][k] = Fy [0][NY-1][k] ;
		dqx[NX][-1][k] = dqx[0][NY-1][k] ;
		dqy[NX][-1][k] = dqy[0][NY-1][k] ;
	}
}
#endif

#if STD
void bounds()
{
	int j ;

	for(j=0;j<NP;j++) {
		/* lower boundary */
		p [-1][j] = p [0][j] ;
		ph [-1][j] = ph [0][j] ;
		F [-1][j] = F [0][j] ;
		dq[-1][j] = dq[0][j] ;

		/* upper boundary */
		p [NN][j] = p [NN-1][j] ;
		ph [NN][j] = ph [NN-1][j] ;
		F [NN][j] = F [NN-1][j] ;
		dq[NN][j] = dq[NN-1][j] ;
	}
}
#endif

