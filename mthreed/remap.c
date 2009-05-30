#include "decs.h"

void remap_linear(REAL (*var)[NY+5][NZ+5], int i, int k, REAL del)
{
        int j,jp,jpp,idel ;
        REAL fshift, tmp[NY] ;

        /* separate into integral & fraction part of shift */
	while(del < 0.) del += NY ;
        idel = (int)floor(del) ;
        fshift = del - idel ;

	/* do integral part of shift */
	for(j=0;j<NY;j++) {
		jp = (j + idel)%NY ; /* gives non-fractional part of (j+idel)/NY */
		jpp = jp + 1 ;
		while(jpp > NY-1) jpp -= NY ;
		tmp[j] = (1. - fshift)*var[i][jp][k] + fshift*var[i][jpp][k] ;
	}
	for(j=0;j<NY;j++) var[i][j][k] = tmp[j] ;
}


void remap_fourier(REAL (*var)[NY+5][NZ+5], int i, int k, REAL del)
{
	REAL tmp[NY] ;
	REAL freq,tmpj,pfac,tfac = 2.*M_PI*del ;
	void drealft(double data[], int n, int isign);
	int j ;
	int Ny = NY ;

	for(j=0;j<NY;j++) tmp[j] = var[i][j][k] ;
	drealft(tmp-1,NY,1) ;

	for(j=2;j<NY;j+=2) {
		freq = ((REAL)j/(2.0*Ny)) ;
		pfac = freq*tfac ;
		tmp[j] = (tmpj=tmp[j])*cos(pfac) - tmp[j+1]*sin(pfac) ;
		tmp[j+1] = tmpj*sin(pfac) + tmp[j+1]*cos(pfac) ;
	}
	pfac = 0.5*tfac ;

	drealft(tmp-1,NY,-1) ;
	for(j=0;j<NY;j++) var[i][j][k] = 2.0*tmp[j]/Ny ;
}


void remap_flux(REAL (*var)[NY+5][NZ+5], int i, int k, int shift, REAL del)
{
	int j,jp,jm,idel,offset ;
	REAL fdel ;
	REAL dq[NY],flux[NY],mdot[NY],trnsvar[NY] ;
	REAL massvar[NY],specvar[NY],fluxvar[NY] ;

        /* separate into integral & fraction part of shift */
        idel = my_nint(del) ;
        fdel = del - (REAL)idel ;

	offset = (var==vy||var==by) ? 1 : 0 ;

	consistent_transport(var,i+shift,k,massvar,specvar,fluxvar,1) ;

	/* calculate mass flux */
	for(j=0;j<NY;j++) trnsvar[j] = fdel ;
	flux_calc(massvar,trnsvar,trnsvar,mdot,offset) ;

	/* calculate flux */
	if(var==rho||var==rhotmp) for(j=0;j<NY;j++) flux[j] = mdot[j] ;
	else flux_calc(specvar,trnsvar,mdot,flux,offset) ;

	/* update transport variable */
	flux_update(fluxvar,flux,offset) ;

	/* do integral part of shift */
        iremap(fluxvar,idel) ;

	/*if(var==vx || var==vy) {
		flux_update(massvar,mdot,offset) ;
		iremap(massvar,idel) ;
	}*/

	consistent_transport(var,i,k,massvar,specvar,fluxvar,-1) ;
}


void iremap(REAL *var, int idel)
{
	static REAL tmp[NY] ;
	int j,jp ;

	for(j=0;j<NY;j++) tmp[j] = var[j] ;

        for(j=0;j<NY;j++) {
                jp = j + idel ;
                while(jp < 0) jp += NY ;
                while(jp > NY-1) jp -= NY ;
                var[jp] = tmp[j] ;
	}
}


void consistent_transport(REAL (*var)[NY+5][NZ+5], int i, int k, REAL *massvar, REAL *specvar, REAL *fluxvar, int fr)
{
	int j,jm ;

	if(fr==1) {
		if(var==vx) for(j=0;j<NY;j++) {
			massvar[j] = 0.5*(rho[i-1][j][k]+rho[i][j][k]) ;
			specvar[j] = var[i][j][k] ;
			fluxvar[j] = massvar[j]*var[i][j][k] ;
		}
		else if(var==vy) for(j=0;j<NY;j++) {
			massvar[j] = 0.5*(rho[i][j-1][k]+rho[i][j][k]) ;
			specvar[j] = var[i][j][k] ;
			fluxvar[j] = massvar[j]*var[i][j][k] ;
		}
		else if(var==bx) for(j=0;j<NY;j++) {
			massvar[j] = 0.5*(rho[i-1][j][k]+rho[i][j][k]) ;
			specvar[j] = var[i][j][k]/massvar[j] ;
			fluxvar[j] = var[i][j][k] ;
		}
		else if(var==by) for(j=0;j<NY;j++) {
			massvar[j] = 0.5*(rho[i][j-1][k]+rho[i][j][k]) ;
			specvar[j] = var[i][j][k]/massvar[j] ;
			fluxvar[j] = var[i][j][k] ;
		}
		else if(var==bz) for(j=0;j<NY;j++) {
			massvar[j] = 0.5*(rho[i][j][k-1]+rho[i][j][k]) ;
			specvar[j] = var[i][j][k]/massvar[j] ;
			fluxvar[j] = var[i][j][k] ;
		}
		else if(var==rho||var==rhotmp) for(j=0;j<NY;j++) {
			massvar[j] = fluxvar[j] = rho[i][j][k] ;
		}
		else for(j=0;j<NY;j++) {
			massvar[j] = rho[i][j][k] ;
			specvar[j] = var[i][j][k]/massvar[j] ;
			fluxvar[j] = var[i][j][k] ;
		}
	}

	else if(fr==-1) {
		if(var==vx) for(j=0;j<NY;j++) var[i][j][k] = 2.0*fluxvar[j]/(rhotmp[i-1][j][k]+rhotmp[i][j][k]) ;
		else if(var==vy) for(j=0;j<NY;j++) {
			jm = j - 1 ;
			while(jm < 0) jm +=NY ;
			var[i][j][k] = 2.0*fluxvar[j]/(rhotmp[i][jm][k]+rhotmp[i][j][k]) ;
		}
  		else for(j=0;j<NY;j++) var[i][j][k] = fluxvar[j] ;
	}
	else {fprintf(stderr,"fr != 1 || -1\n") ; exit(0) ;}
}


void flux_calc(REAL *var, REAL *trnsvar, REAL *mdotvar, REAL *flux, int offset)
{
	int j,jm,jp ;
	REAL dq[NY] ;

	for(j=0;j<NY;j++) {
		jp = j + 1 ;
		while(jp > NY-1) jp -= NY ;
		jm = j - 1 ;
		while(jm < 0)    jm += NY ;
		dq[j] = slope_lim(var[jm],var[j],var[jp]) ;
	}

	for(j=0;j<NY;j++) {
		jp = j + offset ;
		while(jp > NY-1) jp -= NY ;
		jm = jp - 1 ;
		while(jm < 0)    jm += NY ;
		flux[j] = (trnsvar[j] > 0.) ?
			mdotvar[j]*(var[jm] + (1. - trnsvar[j])*dq[jm]) :
			mdotvar[j]*(var[jp] - (1. + trnsvar[j])*dq[jp]) ;
	}
}


void flux_update(REAL *fluxvar, REAL *flux, int offset)
{
	int j,jm,jp ;

	for(j=0;j<NY;j++) {
		jm = j - offset ;
		while(jm < 0)    jm += NY ;
		jp = jm + 1 ;
		while(jp > NY-1) jp -= NY ;
		fluxvar[j] += -(flux[jp]-flux[jm]) ;
	}
}


REAL slope_lim(REAL y1,REAL y2,REAL y3)
{
        REAL Dqm,Dqp,s ;

        /* van leer slope limiter */
	/* N.B.: factor of 2 is missing because factor is
	   absorbed into expressions for fluxes above */
	Dqp = (y3 - y2) ;
	Dqm = (y2 - y1) ;
	s = Dqm*Dqp ;
	if(s <= 0.) return 0. ;
	else return(s/(Dqm+Dqp)) ;

}


int my_nint(REAL z)
{
        if(z > 0) return((int)(z+0.5)) ;
        else return((int)(z-0.5)) ;
}


void bxremap(int i, int k, int idel, int fcase, REAL (*wx)[NX])
{
	REAL tmp[NY],(*bx_dqy)[NY+4][NZ+4] ;
	int j,js,jm,jp ;

#if BVANLEER
	bx_dqy = work1 ;
#endif
	for(j=0;j<NY;j++) {
		js = (j - idel)%NY ;
		while(js < 0) js += NY ;
		while(js > NY - 1) js -= NY ;
		tmp[j] = bx[i][js][k] ;
		switch(fcase) {
		case 1: jm = js - 1 ;
			while(jm < 0) jm += NY ; 
			tmp[j] += wx[0][i]*(bx[i][jm][k] - bx[i][js][k]) ;
#if BVANLEER
			tmp[j] += wx[1][i]*(bx_dqy[i][js][k] - bx_dqy[i][jm][k]) ;
#endif				
			break ;
		case 2: jm = js - 1 ;
			while(jm < 0) jm += NY ;
			tmp[j] += wx[0][i]*(bx[i][jm][k] - bx[i][js][k]) ;
#if BVANLEER
			tmp[j] += wx[1][i]*(bx_dqy[i][js][k] - bx_dqy[i][jm][k]) ;
#endif				
			break ;
		case 3: jp = js + 1 ;
			while(jp > NY - 1) jp -= NY ;
			tmp[j] += wx[0][i]*(bx[i][js][k] - bx[i][jp][k]) ;
#if BVANLEER
			tmp[j] += wx[1][i]*(bx_dqy[i][js][k] - bx_dqy[i][jp][k]) ;
#endif
			break ;
		}
	}
	for(j=0;j<NY;j++) bx[i][j][k] = tmp[j] ;
}


void byremap(int i, int k, int idel, int fcase, REAL (*wy)[NX])
{
	int j,js,jm,jp ;
	REAL tmp[NY],(*bx_dqy)[NY+4][NZ+4],(*bz_dqx)[NY+4][NZ+4],(*bz_dqy)[NY+4][NZ+4] ;
	static REAL w = 1. ;
	
#if BVANLEER
	bx_dqy = work1 ;
	bz_dqx = work2 ;
	bz_dqy = work3 ;
#endif
	for(j=0;j<NY;j++) {
		js = (j - idel)%NY ;
		while(js < 0) js += NY ;
		while(js > NY - 1) js -= NY ;
		jm = js - 1 ;
		while(jm < 0) jm += NY ;
		jp = js + 1 ;
		while(jp > NY - 1) jp -= NY ;
		tmp[j] = wy[13][i]*by[i][js][k]
		       + wy[14][i]*(by[i][jm][k] + by[i][jp][k] - by[i][js][k])
		       + wy[15][i]*(bx[i][jm][k] - bx[i][js][k] + bx[i+1][js][k] - bx[i+1][jm][k])
		       + wy[16][i]*(bz[i][jm][k] - bz[i][js][k] + bz[i][js][k+1] - bz[i][jm][k+1]) ;
		switch(fcase) {
		case 1: 
			tmp[j] += wy[0][i]*bx[i][jm][k] + wy[1][i]*bx[i+1][js][k]
				+ wy[2][i]*(bz[i][jm][k+1] - bz[i][jm][k])
				+ wy[3][i]*(bz[i][js][k+1] - bz[i][js][k]) ;
#if BVANLEER
			tmp[j] += wy[4][i]*bx_dqy[i][jm][k] + wy[5][i]*bx_dqy[i+1][js][k]
				+ wy[6][i]*(bz_dqx[i][jm][k+1] - bz_dqx[i][jm][k])
				+ wy[7][i]*(bz_dqx[i][js][k+1] - bz_dqx[i][js][k])
				+ wy[8][i]*(bz_dqy[i][jm][k+1] - bz_dqy[i][jm][k])
				+ wy[9][i]*(bz_dqy[i][js][k+1] - bz_dqy[i][js][k]) ;
#endif
			break ;
		case 2: tmp[j] += wy[0][i]*bx[i][jm][k] + wy[1][i]*bx[i+1][jm][k]
				+ wy[12][i]*(bz[i][jm][k+1] - bz[i][jm][k]) ;
#if BVANLEER
			tmp[j] += wy[4][i]*bx_dqy[i][jm][k] + wy[5][i]*bx_dqy[i+1][jm][k]
				+ wy[10][i]*(bz_dqx[i][jm][k] - bz_dqx[i][jm][k+1])
				+ wy[11][i]*(bz_dqy[i][jm][k] - bz_dqy[i][jm][k+1]) ;
#endif
			break ;
		case 3: tmp[j] += wy[0][i]*bx[i][js][k] + wy[1][i]*bx[i+1][js][k]
				+ wy[12][i]*(bz[i][js][k+1] - bz[i][js][k]) ;
#if BVANLEER
			tmp[j] += wy[4][i]*(-bx_dqy[i][js][k]) + wy[5][i]*bx_dqy[i+1][js][k]
				+ wy[10][i]*(bz_dqx[i][js][k] - bz_dqx[i][js][k+1])
				+ wy[11][i]*(bz_dqy[i][js][k] - bz_dqy[i][js][k+1]) ;
#endif
			break ;
		}
	}
	for(j=0;j<NY;j++) by[i][j][k] = tmp[j] ;
}


void bzremap(int i, int k, int idel, int fcase, REAL (*wz)[NX])
{
	int j,js,jm,jp ;
	REAL tmp[NY],(*bz_dqx)[NY+4][NZ+4],(*bz_dqy)[NY+4][NZ+4] ;

#if BVANLEER
	bz_dqx = work2 ;
	bz_dqy = work3 ;
#endif

	for(j=0;j<NY;j++) {
		js = (j - idel)%NY ;
		while(js < 0) js += NY ;
		while(js > NY - 1) js -= NY ;
		tmp[j] = bz[i][js][k] ;
		switch(fcase) {
		case 1: jp = js + 1 ;
			while(jp > NY - 1) jp -= NY ;
			jm = js - 1 ;
			while(jm < 0) jm += NY ;
			tmp[j] += wz[0][i]*(bz[i][jm][k] - bz[i][js][k])
			        + wz[1][i]*(bz[i][jp][k] - bz[i][js][k]) ;
#if BVANLEER
			tmp[j] += wz[2][i]*(bz_dqx[i][jm][k] - bz_dqx[i][js][k])
				+ wz[3][i]*(bz_dqx[i][jp][k] - bz_dqx[i][js][k])
				+ wz[4][i]*(bz_dqy[i][jm][k] - bz_dqy[i][js][k])
				+ wz[5][i]*(bz_dqy[i][jp][k] - bz_dqy[i][js][k]) ;
#endif
			break ;
		case 2:	jm = js - 1 ;
			while(jm < 0) jm += NY ;
			tmp[j] += wz[8][i]*(bz[i][jm][k] - bz[i][js][k]) ;
#if BVANLEER
			tmp[j] += wz[6][i]*(bz_dqx[i][js][k] - bz_dqx[i][jm][k])
				+ wz[7][i]*(bz_dqy[i][js][k] - bz_dqy[i][jm][k]) ;
#endif
			break ;
		case 3:	jp = js + 1 ;
			while(jp > NY - 1) jp -= NY ;
			tmp[j] += wz[8][i]*(bz[i][js][k] - bz[i][jp][k]) ;
#if BVANLEER
			tmp[j] += wz[6][i]*(bz_dqx[i][jp][k] - bz_dqx[i][js][k])
				+ wz[7][i]*(bz_dqy[i][jp][k] - bz_dqy[i][js][k]) ;
#endif
			break ;
		}
	}
	for(j=0;j<NY;j++) bz[i][j][k] = tmp[j] ;
}


void weights(int i, int fcase, REAL f, REAL sh, REAL (*wx)[NX], REAL (*wy)[NX], REAL (*wz)[NX])
{
	REAL nm,np,cnm,cnp ;
	int n ;

	nm = sh/2. + f ;
	np = sh/2. - f ;
	cnm = nm*nm/sh ;
	cnp = np*np/sh ;

	switch(fcase) {
		case 1: wx[1][i] = nm*(nm - 1.) ;
			wy[5][i] = np*(1. - np)*dy/dx ;
			break ;
		case 2: wx[1][i] = nm*(nm - 1.) ;
			wy[5][i] = -np*(1. + np)*dy/dx ;
			wz[7][i] = sh*sh/12. - f*(f - 1.) ;
			break ;
		case 3: wx[1][i] = -nm*(nm + 1.) ;
			wy[5][i] = np*(1. - np)*dy/dx ;
			wz[7][i] = sh*sh/12. - f*(f + 1.) ;
			break ;
	}
	
	wx[0][i] = nm ;
	
	wz[0][i] = cnm/2. ;
	wz[1][i] = cnp/2. ;
	wz[2][i] = cnm*(f - sh)/(3.*sh) ;
	wz[3][i] = cnp*(f + sh)/(3.*sh) ;
	wz[4][i] = cnm*(1.5 - nm)/3. ;
	wz[5][i] = cnp*(1.5 - np)/3. ;
	wz[6][i] = sh/6. ;
	wz[8][i] = f ;
	
	wy[0][i] = -nm*dy/dx ;
	wy[1][i] = -np*dy/dx ;
	wy[2][i] = wz[0][i]*dy/dz ;
	wy[3][i] = -wz[1][i]*dy/dz ;
	wy[4][i] = wx[1][i]*dy/dx ;
	wy[6][i] = wz[2][i]*dy/dz ;
	wy[7][i] = -wz[3][i]*dy/dz ;
	wy[8][i] = wz[4][i]*dy/dz ;
	wy[9][i] = -wz[5][i]*dy/dz ;
	wy[10][i] = wz[6][i]*dy/dz ;
	wy[11][i] = wz[7][i]*dy/dz ;
	wy[12][i] = f*dy/dz ;
	wy[13][i] = 0.5 ;
	wy[14][i] = 1. - wy[13][i] ;
	wy[15][i] = wy[14][i]*dy/dx ;
	wy[16][i] = wy[14][i]*dy/dz;
}
