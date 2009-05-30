
#include "decs.h"

/* diagnostics subroutine */


void diag(call_code)
int call_code ;
{
	char dfnam[10],ifnam[10] ;
	double ek,eki,eb,ebi,eth,ethi,etot ;
	double e_fin,amp_fin ;
	double eg,egi ;
	double m_fin ;
	double vxa,vya,vza ;
	double cmode_amp,smode_amp ;
	double Bxa,Bya,Bza ;
	double nj ;
	FILE *dump_file ;
	FILE *im_file ;
	int i,j ;
	void dump() ;

	static double e_init ;
	static double m_init ;
	static double amp_init ;
	static double tdump,timage ;
	static FILE *ener_file ;
	static int dump_cnt,im_cnt ;

	if(call_code==0) {
		ener_file = fopen("ener.out","w") ;
		if(ener_file==NULL) {
			fprintf(stderr,"error opening energy output file\n") ;
			exit(1) ;
		}
		tdump = t - 1.e-12 ;
		timage = t - 1.e-12 ;
		dump_cnt = 0 ;
		im_cnt = 0 ;
	}

	/* calculate energies */
	eb = 0. ;
	ek = 0. ;
	eth = 0. ;
	eg = 0. ;
	etot = 0. ;
	cmode_amp = 0. ;
	smode_amp = 0. ;
	for(i=0;i<NX;i++) 
	for(j=0;j<NY;j++) {
		vxa = 0.5*(vx[i+1][j]+vx[i][j]) ;
		vya = 0.5*(vy[i][j+1]+vy[i][j]) ;
		vza = vz[i][j] ;
		Bxa = 0.5*(Bx[i+1][j]+Bx[i][j]) ;
		Bya = 0.5*(By[i][j+1]+By[i][j]) ;
		Bza = Bz[i][j] ;
		eki = 0.5*rho[i][j]*(vxa*vxa + vya*vya + vza*vza) ;
		ebi = 0.5*(Bxa*Bxa + Bya*Bya + Bza*Bza) ;
		egi = 0.5*rho[i][j]*pot[i][j] ;
		if(fabs(gam - 1.) > 1.e-6) ethi = e[i][j] ;
		else ethi = cs*cs*rho[i][j] ;
		ek += eki*a*dx*a*dy ;
		eb += ebi*a*dx*a*dy ;
		eth += ethi*a*dx*a*dy ;
		eg += egi*a*dx*a*dy ;
		cmode_amp += cos(2.*M_PI*(i+0.5)/NX)*rho[i][j]*a*dx*a*dy ;
		smode_amp += sin(2.*M_PI*(i+0.5)/NX)*rho[i][j]*a*dx*a*dy ;
	}
	etot = ek+eth+eb+eg ;

	if(call_code == 0) {
		e_init = etot ;
		m_init = 0. ;
		amp_init = sqrt(cmode_amp*cmode_amp + smode_amp*smode_amp) ;
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) 
			m_init += rho[i][j]*a*dx*a*dy ;
	}
	else if (call_code == 2) {
		e_fin = etot ;
		amp_fin = sqrt(cmode_amp*cmode_amp + smode_amp*smode_amp) ;
		m_fin = 0. ;
		for(i=0;i<NX;i++) 
		for(j=0;j<NY;j++) 
			m_fin += rho[i][j]*a*dx*a*dy ;
		fprintf(stderr,"\n\nEnergy: ini,fin,del: %g %g %g\n",
			e_init,e_fin,(e_fin-e_init)/e_init) ;
		fprintf(stderr,"mass: ini,fin,del: %g %g %g\n",
			m_init,m_fin,(m_fin-m_init)/m_init) ;
		fprintf(stderr,"amp: ini,fin,del: %g %g %g\n",
			amp_init,amp_fin,(amp_fin-amp_init)/amp_init) ;
		fprintf(stderr,"rho: %g\n",rho[0]) ;
	}
	else {
		fprintf(ener_file,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n",
			t,ek,eth,eb,eg,(etot-e_init)/e_init,cmode_amp,smode_amp) ;
		fflush(ener_file) ;
		nj = sqrt(G/M_PI) ;

		/*
		if(eg < -(M_PI*M_PI*nj*nj/6. - 1.) && G > GMIN) {
			fprintf(stderr,"stop by energy criterion\n") ;
			exit(1) ;
		}
		*/
	}


	/* dump at regular intervals */
	if(t >= tdump || call_code == 2) {
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
