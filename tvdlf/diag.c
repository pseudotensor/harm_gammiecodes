
#include "decs.h"

/* diagnostics subroutine */

void diag(call_code)
int call_code ;
{
	char dfnam[10] ;
	int i,j,k ;
	double rho,v,u,ek,eth,etot,e_fin,m_fin ;
	FILE *dump_file ;
	void dump() ;
	double U[NP],px,py,pz,e,rmed,divb,divbmax ;
	void primtoU(double *p, double *U) ;

	static int dump_cnt ;
	static double e_init,m_init,tdump ;
	static FILE *ener_file ;

	static double timage ;
	char ifnam[50] ;
	static int im_cnt ;
	FILE *im_file ;

	if(call_code==0) {
		ener_file = fopen("ener.out","w") ;
		if(ener_file==NULL) {
			fprintf(stderr,"error opening energy output file\n") ;
			exit(1) ;
		}
		tdump = t + SMALL ;
		timage = t + SMALL ;
		dump_cnt = 0 ;
		im_cnt = 0 ;
	}

	/* calculate energies */
	px = py = pz = e = rmed = divbmax = 0. ;
	LOOP {
		primtoU(p[i][j],U) ;
		rmed += U[0]*dV ;
		px += U[UX]*dV ;
		py += U[UY]*dV ;
		pz += U[UZ]*dV ;
		e += U[UU]*dV ;
		/* flux-ct defn */
		divb = fabs( 0.5*(p[i][j][BX] + p[i][j-1][BX]
			- p[i-1][j][BX] - p[i-1][j-1][BX])/dx +
			0.5*(p[i][j][BY] + p[i-1][j][BY]
			- p[i][j-1][BY] - p[i-1][j-1][BY])/dy) ;
		if(divb > divbmax) divbmax = divb ;
	}
	fprintf(stderr,"divbmax: %g\n",divbmax) ;

	if(call_code == 0) {
		e_init = e ;
		m_init = rmed ;
	}
	else if (call_code == 2) {
		e_fin = e ;
		m_fin = rmed ;
		fprintf(stderr,"\n\nEnergy: ini,fin,del: %g %g %g\n",
			e_init,e_fin,(e_fin-e_init)/e_init) ;
		fprintf(stderr,"mass: ini,fin,del: %g %g %g\n",
			m_init,m_fin,(m_fin-m_init)/m_init) ;
	}
	else {
		fprintf(ener_file,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n",
			t,rmed,px,py,pz,e) ;
		fflush(ener_file) ;
	}


	/* dump at regular intervals */
	if(t >= tdump || call_code == 2 || call_code == 0) {
		fprintf(stderr,"dumping!") ;
		/* make regular dump file */
		sprintf(dfnam,"dump%03d",dump_cnt) ;
		dump_file = fopen(dfnam,"w") ;

		if(dump_file==NULL) {
			fprintf(stderr,"error opening dump file\n") ;
			exit(2) ;
		}

		dump(dump_file) ;
		fclose(dump_file) ;

		dump_cnt++ ;
		tdump += DTd ;
	}


        /* make image of density at regular intervals */
        if(t >= timage || call_code == 2) {
                sprintf(ifnam,"im%04d",im_cnt) ;
                im_file = fopen(ifnam,"w") ;
 
                if(im_file==NULL) {
                        fprintf(stderr,"error opening image file\n") ;
                        exit(2) ;
                }
 
                image(im_file) ;
                fclose(im_file) ;
 
                im_cnt++ ;
                timage += DTi ;
        }


}
