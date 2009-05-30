
REAL rhoas[NX+5][NY+5][NZ+5] ;
REAL rhotmpas[NX+5][NY+5][NZ+5] ;
REAL vxas[NX+5][NY+5][NZ+5] ;
REAL vyas[NX+5][NY+5][NZ+5] ;
REAL vzas[NX+5][NY+5][NZ+5] ;
REAL bxas[NX+5][NY+5][NZ+5] ;
REAL byas[NX+5][NY+5][NZ+5] ;
REAL bzas[NX+5][NY+5][NZ+5] ;
REAL eas[NX+5][NY+5][NZ+5] ;
REAL work1as[NX+4][NY+4][NZ+4] ;
REAL work2as[NX+4][NY+4][NZ+4] ;
REAL work3as[NX+4][NY+4][NZ+4] ;
REAL work4as[NX+4][NY+4][NZ+4] ;
REAL work5as[NX+4][NY+4][NZ+4] ;
REAL work6as[NX+4][NY+4][NZ+4] ;

REAL (*rho)[NY+5][NZ+5] ;
REAL (*rhotmp)[NY+5][NZ+5] ;
REAL (*vx)[NY+5][NZ+5] ;
REAL (*vy)[NY+5][NZ+5] ;
REAL (*vz)[NY+5][NZ+5] ;
REAL (*bx)[NY+5][NZ+5] ;
REAL (*by)[NY+5][NZ+5] ;
REAL (*bz)[NY+5][NZ+5] ;
REAL (*e)[NY+5][NZ+5] ;
REAL (*work1)[NY+4][NZ+4] ;
REAL (*work2)[NY+4][NZ+4] ;
REAL (*work3)[NY+4][NZ+4] ;
REAL (*work4)[NY+4][NZ+4] ;
REAL (*work5)[NY+4][NZ+4] ;
REAL (*work6)[NY+4][NZ+4] ;

REAL cour ;
REAL cs ;
REAL c ;
REAL gam ;
REAL dx,dy,dz ;
REAL dt ;
REAL t,tf ;
REAL nu_vnr ;
REAL nu_l ;

REAL Lx ;
REAL Ly ;
REAL Lz ;

REAL dV ;

REAL DTd ;
REAL DTi ;
REAL DTl ;

REAL drive ;
REAL A0,B0,B0x,B0y,B0z ;
REAL work,workp ;
REAL W,Wz,q ;

REAL kx,ky,kz ;
REAL remr[8],remi[8],norm,pert ;

int nstep ;

