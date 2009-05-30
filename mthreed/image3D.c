/* 
	produces an "r8" file.  
*/

#include "decs.h"

#define DIVB(i,j,k) ((bx[(i)+1][(j)][(k)]-bx[(i)][(j)][(k)])/dx+(by[(i)][(j)+1][(k)]-by[(i)][(j)][(k)])/dy+(bz[(i)][(j)][(k)+1]-bz[(i)][(j)][(k)])/dz)

void image(FILE *fp)
{
	REAL iq,liq,a,b,lmax,lmin ;
	int i,j ;
	int size ;

	if(NX <= NY && NX <= NZ) size = NX ;
	else if(NY <= NX && NY <= NZ) size = NY ;
	else size = NZ ;

	/* density mapping is logarithmic, in 255 steps
	   between e^lmax and e^lmin */

	lmax = 2.5 ;
	lmin = -0.5 ;
#ifdef LINEAR
	lmax = pert ;
	lmin = -pert ;
#endif
	a = 256./(lmax - lmin) ;
	b = -a*lmin ;
	
	/* N.B.: this only really works for cubes */
	for(j=(NY + size)/2-1;j>= (NY - size)/2;j--) {
		for(i=(NX-size)/2;i<(NX+size)/2;i++) {
			iq = DIVB(i,j,NZ/2)/4. ;
			iq = vy[i][j][NZ/2]  ;
			liq = (iq/cs)*128 + 128 ;
			if(liq > 255.) liq = 255. ;
			if(liq < 0.) liq = 0. ;
			fprintf(fp,"%c",(char)((int)liq)) ;
		}
		fprintf(fp,"%c",0) ;
		fprintf(fp,"%c",0) ;
		for(i=(NX-size)/2;i<(NX+size)/2;i++) {
			iq = rho[i][j][NZ/2] ;
			liq = ((iq - 1.)/1.)*128 + 128 ;
			if(liq > 255.) liq = 255. ;
			if(liq < 0.) liq = 0. ;
			fprintf(fp,"%c",(char)((int)liq)) ;
		}
	}
	for(i=0;i<2*size+2;i++) fprintf(fp,"%c",0) ;
	for(i=0;i<2*size+2;i++) fprintf(fp,"%c",0) ;
	for(j=(NZ + size)/2-1;j>= (NZ - size)/2;j--) {
		for(i=(NX-size)/2;i<(NX+size)/2;i++) {
			iq = DIVB(i,NY/2,j)/4. ;
			iq = vy[i][NY/2][j] ;
			liq = (iq/cs)*128 + 128 ;
			if(liq > 255.) liq = 255. ;
			if(liq < 0.) liq = 0. ;
			fprintf(fp,"%c",(char)((int)liq)) ;
		}
		fprintf(fp,"%c",0) ;
		fprintf(fp,"%c",0) ;
		for(i=(NX-size)/2;i<(NX+size)/2;i++) {
			iq = rho[i][NY/2][j] ;
			liq = ((iq - 1.)/1.)*128 + 128 ;
			if(liq > 255.) liq = 255. ;
			if(liq < 0.) liq = 0. ;
			fprintf(fp,"%c",(char)((int)liq)) ;
		}
	}
	for(i=0;i<2*size+2;i++) fprintf(fp,"%c",0) ;
	for(i=0;i<2*size+2;i++) fprintf(fp,"%c",0) ;
	for(j=(NZ + size)/2-1;j>= (NZ - size)/2;j--) {
		for(i=(NY-size)/2;i<(NY+size)/2;i++) {
			iq = DIVB(NX/2,i,j)/4. ;
			iq = vy[NX/2][i][j] ;
			liq = (iq/cs)*128 + 128 ;
			if(liq > 255.) liq = 255. ;
			if(liq < 0.) liq = 0. ;
			fprintf(fp,"%c",(char)((int)liq)) ;
		}
		fprintf(fp,"%c",0) ;
		fprintf(fp,"%c",0) ;
		for(i=(NY-size)/2;i<(NY+size)/2;i++) {
			iq = rho[NX/2][i][j] ;
			liq = ((iq - 1.)/1.)*128 + 128 ;
			if(liq > 255.) liq = 255. ;
			if(liq < 0.) liq = 0. ;
			fprintf(fp,"%c",(char)((int)liq)) ;
		}
	}
}
