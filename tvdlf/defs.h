

/* primitive variables */
double a_p[NX+2][NY+2][NP] ;
double a_dqx[NX+2][NY+2][NP] ;
double a_dqy[NX+2][NY+2][NP] ;
double a_Fx[NX+2][NY+2][NP] ;
double a_Fy[NX+2][NY+2][NP] ;
double a_ph[NX+2][NY+2][NP] ;

double (*  p)[NY+2][NP] ;
double (*dqx)[NY+2][NP] ;
double (*dqy)[NY+2][NP] ;
double (* Fx)[NY+2][NP] ;
double (* Fy)[NY+2][NP] ;
double (* ph)[NY+2][NP] ;
double (*emf)[NY+2][NP] ;

double Lx,Ly ;

double cour ;
double cs ;
double gam ;
double dx,dy,dV ;
double dt ;
double t,tf ;

double DTd ;
double DTl ;
double DTi ;

