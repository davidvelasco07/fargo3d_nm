#include "fargo3d.h"

void ReMapGeometryStock(real inner, real outer) { 

  OUTPUT2D(Vx0);
  OUTPUT2D(Vy0);
  OUTPUT2D(Density0);
  OUTPUT2D(Energy0);

  real* vx0 = Vx0->field_cpu;
  real* vy0 = Vy0->field_cpu;
  real* rho0 = Density0->field_cpu;
  real* e0   = Energy0->field_cpu;
  real r;
  real omega;

  real dy;

  int i,j,k;

  int size_y= Ny+2*NGHY;
  int size_z= Nz+2*NGHZ;

  YMIN = inner;
  YMAX = outer;
  
  dy = (YMAX-YMIN)/NY;
  
  for (j=0; j<Ny+2*NGHY+1; j++) {    
    Ymin(j) = YMIN + dy*(j+Y0-NGHY);
  }
  for (j = 0; j<Ny+2*NGHY; j++) {
    Ymed(j) = 0.5*(Ymin(j+1)+Ymin(j));
  }
  for (j = 1; j<Ny+2*NGHY; j++) {
    InvDiffYmed(j) = 1./(Ymed(j)-Ymed(j-1));
  }

  //Updating Stockholm

  j = k = 0;

#ifdef Z
  for (k=0; j<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
    r = Ymed(j);
    omega = sqrt(G*MSTAR/r/r/r);    
    rho0[l2D] = SIGMA0*pow((r/R0),-SIGMASLOPE);
    e0[l2D]   = sqrt(G*MSTAR/r)*ASPECTRATIO*pow(r/R0, FLARINGINDEX);
    vx0[l2D]  = omega*r*sqrt(1.0+pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*(2.0*FLARINGINDEX - 1.0 - SIGMASLOPE)) - OMEGAFRAME*r;
    vy0[l2D]  = 0.0;
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
  
  InitSurfaces();

#ifdef GPU
  //Updating necessary lenghts, surfaces and volumes in the GPU
  DevMemcpyH2D(Ymin_d,Ymin,sizeof(real)*(Ny+2*NGHY+1));
  DevMemcpyH2D(Sxj_d,Sxj,sizeof(real)*(Ny+2*NGHY));
  DevMemcpyH2D(Syj_d,Syj,sizeof(real)*(Ny+2*NGHY));
  DevMemcpyH2D(Szj_d,Szj,sizeof(real)*(Ny+2*NGHY));
  DevMemcpyH2D(InvVj_d,InvVj,sizeof(real)*(Ny+2*NGHY));
#endif

}
