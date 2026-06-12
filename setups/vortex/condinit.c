#include "fargo3d.h"

void InitDensEnergy() {

  int i,j,k;
  real *rho, *energy;

  rho = Density->field_cpu;
  energy = Energy->field_cpu;
    
  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = 0;
  int end_i = Nx;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
	rho[l] = 1.0;
	energy[l] = 1.0/(GAMMA*(GAMMA-1.0)*MACHNUMBER*MACHNUMBER);
#ifdef ISOTHERMAL
	energy[l] = 1.0/MACHNUMBER;
#endif
      }
    }
  }
}

void InitVazimPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t;
  real rho;
  real vt, omega;
  real *vrad;
  real *vphi;

  vphi = Vx->field_cpu;
  vrad = Vy->field_cpu;
    
  boolean GhostInclude = TRUE;

  real x,y,h,factor,vr;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = 0;
  int end_i = Nx;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
	r = Ymed(j);
	omega = sqrt(G*MSTAR/r/r/r);
	vt = omega*r;
	vt -= OMEGAFRAME*r;

	/* Vphi centering below */
	x = Ymed[j]*cos(Xmin[i])-sqrt(.5);
	y = Ymed[j]*sin(Xmin[i])-sqrt(.5);
	h = 0.1/2.;
	factor = -exp(-(x*x+y*y)/(h*h));
	vt += factor*(sin(Xmin[i])*y+cos(Xmin[i])*x);

	/* Vrad centering below */
	x = Ymin[j]*cos(Xmed[i])-sqrt(.5);
	y = Ymin[j]*sin(Xmed[i])-sqrt(.5);
	h = 0.1/2.;
	factor = -exp(-(x*x+y*y)/(h*h));
	vr = factor*(cos(Xmed[i])*(-y)+sin(Xmed[i])*x);

	vphi[l] = vt;
	vrad[l] = vr;
      }
    }
  }    
}

void CondInit() {
  
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);

  InitDensEnergy ();
  InitVazimPlanet ();
}
