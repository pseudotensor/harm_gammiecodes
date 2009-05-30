
#include "decs.h"

double Uc[NP] ;
double Ut[NP-3] ;
double Bt[3] ;

#define NMAX	200

/* pr *MUST* contain initial guess */
void Utoprim(double *U, double *pr)
{
	double tolx,tolf,pr_save[NP],U0[NP] ;
	int N,ntrial,i,k,test ;
	int mnewt(int ntrail, double *p, int n, double tolx, double tolf) ;
	void primtoU(double *pr, double *U) ;
	void diag(int call_code) ;


	/* polish up the guess w/ Newton-Raphson */
	ntrial = 30 ;
	tolx = 1.e-11;
	tolf = 1.e-11 ;
	N = 1 ;
	primtoU(pr,U0) ;
	PLOOP pr_save[k] = pr[k] ;
	while( 1 ) {
		PLOOP pr[k] = pr_save[k] ;
		for(i=1;i<=N;i++) {
			PLOOP Uc[k] = ((N-i)*U0[k] + i*U[k])/N ;
			Bpare(Uc,Bt) ;
			tpare(U0,guess) ;
			tpare(Uc,Ut) ;
			test = mnewt(ntrial,guess-1,NP-3, tolx, tolf) ;
			if( !test ) break ;
		}
                if(test) {
			join(pr,guess,Bt) ;
                        return ;
                }
                else N++ ;
                fprintf(stderr,"N,i: %d %d\n",N,i) ;
                if(N > NMAX) {
                        diag(2) ;
                        exit(1) ;
                }
	}
}

void Bpare(double *U,*Bo)
{
	Bo[0] = U[4] ;
	Bo[1] = U[5] ;
	Bo[2] = U[6] ;
}
void tpare(double *U,double *to,double *Bo)
{
	to[0] = U[0] ;
	to[1] = U[1] ;
	to[2] = U[2] ;
	to[3] = U[3] ;
	to[4] = U[7] ;
}
void join(double *pn,double *to,double *Bo)
{
	pn[0] = to[0] ;
	pn[1] = to[1] ;
	pn[2] = to[2] ;
	pn[3] = to[3] ;
	pn[4] = Bo[0] ;
	pn[5] = Bo[1] ;
	pn[6] = Bo[2] ;
	pn[7] = to[4] ;
}

/* auxiliary function required by mnewt */
void usrfun(double *tr,int n,double *beta,double **alpha)
{
	double Ut[NP] ;
	int k ;
	void primtoU(double *p, double *U), dudp_calc(double *p, double **a),
		ndudp_calc(double *p, double **a);

	join(pr,tr+1,Bt) ;
	primtoU(pr,Ut) ;
	beta[1] = Ut[0] - tr[0] ;
	beta[2] = Ut[1] - tr[1] ;
	beta[3] = Ut[2] - tr[2] ;
	beta[4] = Ut[3] - tr[3] ;
	beta[5] = Ut[7] - tr[4] ;
	dudp_calc(tr+1,alpha) ;
}

void dudp_calc(double *pr, double **dudp) 
{
	double r,ux,uy,uz,Bx,By,Bz,P,U,w,u0 ;
	double b0,bx,by,bz,b2 ;
	double db0dBx,db0dBy,db0dBz,db0dux,db0duy,db0duz ;
	double dbxdBx,dbxdBy,dbxdBz,dbydBx,dbydBy,dbydBz,dbzdBx,dbzdBy,dbzdBz ;
	double db2dBx,db2dBy,db2dBz,du0dux,du0duy,du0duz ;
	double dbxdux,dbxduy,dbxduz,dbydux,dbyduy,dbyduz,dbzdux,dbzduy,dbzduz ;
	double db2dux,db2duy,db2duz ;
	int j,k ;

	r = pr[RHO] ;
	ux = pr[UX] ;
	uy = pr[UY] ;
	uz = pr[UZ] ;
	Bx = pr[BX] ;
	By = pr[BY] ;
	Bz = pr[BZ] ;
	U = pr[UU] ;

	P = (gam - 1.)*U ;
	w = P + r + U ;
	u0 = sqrt(1. + ux*ux + uy*uy + uz*uz) ;

	b0 = Bx*ux + By*uy + Bz*uz ;
	bx = (Bx + b0*ux)/u0 ;
	by = (By + b0*uy)/u0 ;
	bz = (Bz + b0*uz)/u0 ;
	b2 = -b0*b0 + bx*bx + by*by + bz*bz ;

	db0dBx = ux ;
	db0dBy = uy ;
	db0dBz = uz ;

	db0dux = Bx ;
	db0duy = By ;
	db0duz = Bz ;

	dbxdBx = (1. + ux*db0dBx)/u0 ;
	dbxdBy = ux*db0dBy/u0 ;
	dbxdBz = ux*db0dBz/u0 ;

	dbydBx = uy*db0dBx/u0 ;
	dbydBy = (1. + uy*db0dBy)/u0 ;
	dbydBz = uy*db0dBz/u0 ;

	dbzdBx = uz*db0dBx/u0 ;
	dbzdBy = uz*db0dBy/u0 ;
	dbzdBz = (1. + uz*db0dBz)/u0 ;

	db2dBx = -2.*b0*db0dBx + 2.*bx*dbxdBx + 2.*by*dbydBx + 2.*bz*dbzdBx ;
	db2dBy = -2.*b0*db0dBy + 2.*bx*dbxdBy + 2.*by*dbydBy + 2.*bz*dbzdBy ;
	db2dBz = -2.*b0*db0dBz + 2.*bx*dbxdBz + 2.*by*dbydBz + 2.*bz*dbzdBz ;

	du0dux = ux/u0 ;
	du0duy = uy/u0 ;
	du0duz = uz/u0 ;

	dbxdux = (b0 + ux*db0dux)/u0 - bx*du0dux/u0 ;
	dbxduy = (ux*db0duy)/u0 - bx*du0duy/u0 ;
	dbxduz = (ux*db0duz)/u0 - bx*du0duz/u0 ;

	dbydux = (uy*db0dux)/u0 - by*du0dux/u0 ;
	dbyduy = (b0 + uy*db0duy)/u0 - by*du0duy/u0 ;
	dbyduz = (uy*db0duz)/u0 - by*du0duz/u0 ;

	dbzdux = (uz*db0dux)/u0 - bz*du0dux/u0 ;
	dbzduy = (uz*db0duy)/u0 - bz*du0duy/u0 ;
	dbzduz = (b0 + uz*db0duz)/u0 - bz*du0duz/u0 ;

	db2dux = -2.*b0*db0dux + 2.*bx*dbxdux + 2.*by*dbydux + 2.*bz*dbzdux ; 
	db2duy = -2.*b0*db0duy + 2.*bx*dbxduy + 2.*by*dbyduy + 2.*bz*dbzduy ; 
	db2duz = -2.*b0*db0duz + 2.*bx*dbxduz + 2.*by*dbyduz + 2.*bz*dbzduz ; 

	dudp[1][1] = u0 ;
	dudp[1][2] = r*du0dux ;
	dudp[1][3] = r*du0duy ;
	dudp[1][4] = r*du0duz ;
	dudp[1][5] = 0. ;
	dudp[1][6] = 0. ;
	dudp[1][7] = 0. ;
	dudp[1][8] = 0. ;

	dudp[2][1] = u0*ux ;
	dudp[2][2] = db2dux*u0*ux + (w + b2)*du0dux*ux + (w + b2)*u0 - db0dux*bx - b0*dbxdux ;
	dudp[2][3] = db2duy*u0*ux + (w + b2)*du0duy*ux  - db0duy*bx - b0*dbxduy ;
	dudp[2][4] = db2duz*u0*ux + (w + b2)*du0duz*ux  - db0duz*bx - b0*dbxduz ;
	dudp[2][5] = db2dBx*u0*ux - db0dBx*bx - b0*dbxdBx ;
	dudp[2][6] = db2dBy*u0*ux - db0dBy*bx - b0*dbxdBy ;
	dudp[2][7] = db2dBz*u0*ux - db0dBz*bx - b0*dbxdBz ;
	dudp[2][8] = gam*u0*ux ;


	dudp[3][1] = u0*uy ;
	dudp[3][2] = db2dux*u0*uy + (w + b2)*du0dux*uy - db0dux*by - b0*dbydux ;
	dudp[3][3] = db2duy*u0*uy + (w + b2)*du0duy*uy + (w + b2)*u0 - db0duy*by - b0*dbyduy ;
	dudp[3][4] = db2duz*u0*uy + (w + b2)*du0duz*uy  - db0duz*by - b0*dbyduz ;
	dudp[3][5] = db2dBx*u0*uy - db0dBx*by - b0*dbydBx ;
	dudp[3][6] = db2dBy*u0*uy - db0dBy*by - b0*dbydBy ;
	dudp[3][7] = db2dBz*u0*uy - db0dBz*by - b0*dbydBz ;
	dudp[3][8] = gam*u0*uy ;

	dudp[4][1] = u0*uz ;
	dudp[4][2] = db2dux*u0*uz + (w + b2)*du0dux*uz - db0dux*bz - b0*dbzdux ;
	dudp[4][3] = db2duy*u0*uz + (w + b2)*du0duy*uz - db0duy*bz - b0*dbzduy ;
	dudp[4][4] = db2duz*u0*uz + (w + b2)*du0duz*uz + (w + b2)*u0 - db0duz*bz - b0*dbzduz ;
	dudp[4][5] = db2dBx*u0*uz - db0dBx*bz - b0*dbzdBx ;
	dudp[4][6] = db2dBy*u0*uz - db0dBy*bz - b0*dbzdBy ;
	dudp[4][7] = db2dBz*u0*uz - db0dBz*bz - b0*dbzdBz ;
	dudp[4][8] = gam*u0*uz ;

	dudp[5][1] = 0. ;
	dudp[5][2] = 0. ;
	dudp[5][3] = 0. ;
	dudp[5][4] = 0. ;
	dudp[5][5] = 1. ;
	dudp[5][6] = 0. ;
	dudp[5][7] = 0. ;
	dudp[5][8] = 0. ;

	dudp[6][1] = 0. ;
	dudp[6][2] = 0. ;
	dudp[6][3] = 0. ;
	dudp[6][4] = 0. ;
	dudp[6][5] = 0. ;
	dudp[6][6] = 1. ;
	dudp[6][7] = 0. ;
	dudp[6][8] = 0. ;

	dudp[7][1] = 0. ;
	dudp[7][2] = 0. ;
	dudp[7][3] = 0. ;
	dudp[7][4] = 0. ;
	dudp[7][5] = 0. ;
	dudp[7][6] = 0. ;
	dudp[7][7] = 1. ;
	dudp[7][8] = 0. ;

	dudp[8][1] = u0*(u0 - 1.) ;
	dudp[8][2] = du0dux*((w + b2)*u0 - r) + u0*(db2dux*u0 + (w + b2)*du0dux) 
		- 0.5*db2dux - 2.*b0*db0dux ;
	dudp[8][3] = du0duy*((w + b2)*u0 - r) + u0*(db2duy*u0 + (w + b2)*du0duy) 
		- 0.5*db2duy - 2.*b0*db0duy ;
	dudp[8][4] = du0duz*((w + b2)*u0 - r) + u0*(db2duz*u0 + (w + b2)*du0duz) 
		- 0.5*db2duz - 2.*b0*db0duz ;
	dudp[8][5] = db2dBx*u0*u0 - 0.5*db2dBx - 2.*b0*db0dBx ;
	dudp[8][6] = db2dBy*u0*u0 - 0.5*db2dBy - 2.*b0*db0dBy ;
	dudp[8][7] = db2dBz*u0*u0 - 0.5*db2dBz - 2.*b0*db0dBz ;
	dudp[8][8] = gam*u0*u0 - (gam - 1.) ;

}

/* this is not covariant */
void field_cd(double *po, double *pn, double *emfp) 
{
	double uxo,uyo,uzo,uxn,uyn,uzn,u0o,u0n,Bxo,Byo,Bzo,Bxn,Byn,Bzn ;
	double vxo,vyo,vzo,vxn,vyn,vzn ;

	uxo = po[UX] ;
	uyo = po[UY] ;
	uzo = po[UZ] ;
	u0o = sqrt(1. + uxo*uxo + uyo*uyo + uzo*uzo) ;

	uxn = pn[UX] ;
	uyn = pn[UY] ;
	uzn = pn[UZ] ;
	u0n = sqrt(1. + uxn*uxn + uyn*uyn + uzn*uzn) ;

	vxo = uxo/u0o ;
	vyo = uyo/u0o ;
	vzo = uzo/u0o ;

	vxn = uxn/u0n ;
	vyn = uyn/u0n ;
	vzn = uzn/u0n ;

	Bxo = po[BX] ;
	Byo = po[BY] ;
	Bzo = po[BZ] ;
	Bxn = pn[BX] ;
	Byn = pn[BY] ;
	Bzn = pn[BZ] ;

	/* EMFs */
	emfp[BX] = -0.5*( vyo*Bzo - vzo*Byo + vyn*Bzn - vzn*Byn ) ;
	emfp[BY] = -0.5*( vzo*Bxo - vxo*Bzo + vzn*Bxn - vxn*Bzn ) ;
	emfp[BZ] = -0.5*( vxo*Byo - vyo*Bxo + vxn*Byn - vyn*Bxn ) ;
}
