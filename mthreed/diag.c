
#include "decs.h"

/* diagnostics subroutine */

void diag(call_code)
int call_code ;
{
	char dfnam[10],ifnam[10] ;
	REAL ekx,eky,ekz,eg,eth,wm,wr,divb,divbmax,mass ;
	REAL ebx,eby,ebz,vxbar,vybar ;
	REAL m_fin ;
	REAL vxa,vya,vza,bxa,bya,bza ;
	REAL X,Y,Z ;
	FILE *dump_file ;
	FILE *im_file ;
	int i,j,k ;
	void dump() ;

	static REAL m_init ;
	static REAL tdump,timage ;
	static FILE *ener_file ;
	static int dump_cnt,im_cnt ;

	/* write to logfile */
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
	ekx = 0. ;
	eky = 0. ;
	ekz = 0. ;
	ebx = 0. ;
	eby = 0. ;
	ebz = 0. ;
	eth = 0. ;
	eg = 0. ;
	wm = 0. ;
	wr = 0. ;
	divb = 0. ;

	/* find mean x and y velocities */
	mass = vxbar = vybar = 0. ;
	LOOP {
		mass += rho[i][j][k] ;
		vxa = 0.5*(vx[i+1][j][k]+vx[i][j][k]) ;
		vxbar += rho[i][j][k]*vxa ;
		vya = 0.5*(vy[i][j+1][k]+vy[i][j][k]) ;
		vybar += rho[i][j][k]*vya ;
	}
	vxbar /= mass ;
	vybar /= mass ;

	ekx = eky = ekz = ebx = eby = ebz = eth = eg = 0. ;
	wm = wr = divbmax = 0. ;
	LOOP {
		/* average to center of zone */
		X = (i + 0.5)*dx - 0.5*Lx ;
		Z = (i + 0.5)*dz - 0.5*Lz ;

		vxa = 0.5*(vx[i+1][j][k]+vx[i][j][k]) ;
		vya = 0.5*(vy[i][j+1][k]+vy[i][j][k]) ;
		vza = 0.5*(vz[i][j][k+1]+vz[i][j][k]) ;
		bxa = 0.5*(bx[i+1][j][k]+bx[i][j][k]) ;
		bya = 0.5*(by[i][j+1][k]+by[i][j][k]) ;
		bza = 0.5*(bz[i][j][k+1]+bz[i][j][k]) ;

		ekx += 0.5*rho[i][j][k]*vxa*vxa ;
		eky += 0.5*rho[i][j][k]*vya*vya ;
		ekz += 0.5*rho[i][j][k]*vza*vza ;
		ebx += 0.5*bxa*bxa ;
		eby += 0.5*bya*bya ;
		ebz += 0.5*bza*bza ;

		eth += cs*cs*rho[i][j][k]*log(rho[i][j][k]) ;

		eg += 0.5*rho[i][j][k]*Z*Z*Wz*Wz ;

		wm += -bxa*bya ;
		wr += rho[i][j][k]*vxa*vya ;

		divb = fabs(
			(bx[i+1][j][k] - bx[i][j][k])/dx
		    + (by[i][j+1][k] - by[i][j][k])/dy
		    + (bz[i][j][k+1] - bz[i][j][k])/dz 
		    ) ;
		if(divb > divbmax) divbmax = divb ;
	}
	ekx /= NX*NY*NZ ;
	eky /= NX*NY*NZ ;
	ekz /= NX*NY*NZ ;
	ebx /= NX*NY*NZ ;
	eby /= NX*NY*NZ ;
	ebz /= NX*NY*NZ ;
	eth /= NX*NY*NZ ;
	eg /= NX*NY*NZ ;
	wm /= NX*NY*NZ ;
	wr /= NX*NY*NZ ;

	if(call_code == 0) {
		m_init = mass ;
	}
	else if (call_code == 2) {
		m_fin = mass ;

		fprintf(stderr,"mass: ini,fin,del: %g %g %g\n",
			m_init,m_fin,(m_fin-m_init)/m_init) ;
	}
	else {
		fprintf(ener_file,
			"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g ",
			t,ekx,eky,ekz,ebx,eby,ebz,eth,eg,wm,wr,divbmax) ;
		/* supplemental diagnostics here */
		fprintf(ener_file,"%10.5g %10.5g ",vxbar,vybar) ;

		fprintf(ener_file,"\n") ;

		fflush(ener_file) ;
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
                sprintf(ifnam,"im%04d.ppm",im_cnt) ;
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
