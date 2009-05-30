
/* test brio and wu solution against results of Ryu
   and Jones, Table 5A */

#include "decs.h"

void init()
{
	int i,j,k ;
	void set_arrays() ;
	void bounds() ;
	double X,Y,R,rho0,amp ;
	double frho(double R), fp(double R) ;
	int smooth ;
	FILE * in;
	double ftemp[NP],fnull;

	/* assume throughout that G = rho = L = 1 */
	gam = 5./3. ;
	dt = 0.001 ;
	Lx = 2.*M_PI ;
	Ly = 2.*M_PI ;

	/* set up arrays,details */
	dx = Lx/NX ;
	dy = Ly/NY ;
	dV = dx*dy ;
	t = 0. ;
	set_arrays() ;

	/* more hardwired choices */
	tf = 8.0 ;
	DTd = tf/10. ;	/* dumping frequency */
	DTi = 0.05 ;
	DTl = 0.01 ;	/* logfile frequency */

	amp = 1.e-4 ;
	rho0 = 1. ;

#if 0
	/* Orszag-Tang vortex */
	cour = 0.25 ;

	LOOP {
		X = (i + 0.5)*dx - 0.5*Lx ;
		Y = (j + 0.5)*dy - 0.5*Ly ;

		p[i][j][UX] = -sin(Y + M_PI) ;
		p[i][j][UY] = sin(X + M_PI) ;
		p[i][j][UZ] = 0. ;
		p[i][j][BX] = -sin(Y + M_PI) ;
		p[i][j][BY] = sin(2.*(X + M_PI)) ;
		p[i][j][BZ] = 0. ;

		p[i][j][RHO] = 25./9. ;
		p[i][j][UU] = (5./3.)/(gam - 1.) ;
	}
	
	if( (in=fopen("../dump0019.dat","rt"))==NULL){
	  printf("can't open dump\n");
	  fflush(stdout);
	  exit(1);
	}
	else{
	  while(fgetc(in)!='\n');
	  while(fgetc(in)!='\n');
	  while(fgetc(in)!='\n');
	  while(fgetc(in)!='\n');
	  while(fgetc(in)!='\n');
	  while(fgetc(in)!='\n');
	  while(fgetc(in)!='\n');
	  for(j=-2;j<NY+2;j++){
	    for(i=-2;i<NX+2;i++){
	      fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&ftemp[0],&ftemp[7],&fnull,&ftemp[1],&ftemp[2],&ftemp[3],&ftemp[4],&ftemp[5],&ftemp[6]);
	      while(fgetc(in)!='\n');
	      if((i>=-1)&&(i<=NX)&&(j>=-1)&&(j<=NY) ){
		PLOOP{
		  p[i][j][k]=ftemp[k];
		}
	      }
	    }
	  }
	}

#endif
#if 1
	/* Komissarov's explosion problem */
	tf = 4.0 ;
	cour = 0.4 ;

	LOOP {
		X = (i + 0.5)*dx - 0.5*Lx ;
		Y = (j + 0.5)*dy - 0.5*Ly ;
		R = sqrt(X*X + Y*Y) ;

		p[i][j][UX] = 0. ;
		p[i][j][UY] = 0. ;
		p[i][j][UZ] = 0. ;
		p[i][j][BX] = 2.0 - 1.5*cos(2.*M_PI*Y/Ly) ;
		p[i][j][BY] = 0. ;
		p[i][j][BZ] = 0. ;

		p[i][j][RHO] = frho(R) ;
		p[i][j][UU] = fp(R)/(gam - 1.) ;
	}
#endif

	/* enforce boundary conditions */
	bounds() ;

}

double frho(double R)
{
	double Rout,Rin,dR,rmax,rmin ;

	/*
	rmax = 1.e-2 ;
	rmin = 1.e-4 ;
	*/
	rmax = 1.e-1 ;
	rmin = 1.e-3 ;
	Rout = 0.8 ;
	Rin = 0.8 ;
	dR = (Rout - Rin) ;
	
	if(R > Rout) return(rmin) ;
	if(R < Rin) return(rmax) ;
	else {
		return( exp(
			log(rmax)*(Rout - R)/dR +
			log(rmin)*(R - Rin)/dR)) ;
	}
}
double fp(double R) 
{
	double Rout,Rin,dR,pmin,pmax ;

	/*
	pmin = 3.e-5 ;
	pmax = 1. ;
	*/
	pmin = 1.e-2 ;
	pmax = 1. ;
	Rout = 0.8 ;
	Rin = 0.8 ;
	dR = (Rout - Rin) ;

	if(R > Rout) return(pmin) ;
	else if(R < Rin) return(pmax) ;
	else {
		return( exp(
			log(pmax)*(Rout - R)/dR +
			log(pmin)*(R - Rin)/dR)) ;
	}

}
