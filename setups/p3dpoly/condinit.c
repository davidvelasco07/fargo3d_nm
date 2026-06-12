#include "fargo3d.h"
#include <J_jupiter.h>
void CondInit() {
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Pot);
  OUTPUT(Vx);
  OUTPUT(Vy);
  OUTPUT(Vz);
  int i,j,k;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *rho;
  real h;
  real ztop, rhomidplane,z,vphis;
  real omega,K;
  real r, r3;
  printf("Initializing values for grid %d of level %d\n",Density->desc->number, Density->level);
  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      h = ASPECTRATIO*Ymed(j);
      r = Ymed(j);
      r3 = r*r*r;
      omega = sqrt(G*MSTAR/(r3));
      ztop = sqrt(2.0/(GAMMA-1.0))*ASPECTRATIO*r;
      for (i=NGHX; i<Nx+NGHX; i++) {
	v2[l] = v3[l] = 0.0;
	z = r*cos(Zmed(k));
	rhomidplane = pow(r/R0,-2.5);
	rho[l] = rhomidplane*pow((1.-z*z/ztop/ztop),1./(GAMMA-1.0));
	e[l] = K = ASPECTRATIO*ASPECTRATIO*R0*R0/GAMMA; //OMEGA_K at R0=1
      }
    }
  }
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      r = Ymed(j);
      for (i=NGHX; i<Nx+NGHX; i++) {
	vphis = (G*MSTAR/r);
	if ((j > 0) && (j < Ny+2*NGHY-1))
	  vphis = vphis+(K*pow(rho[lyp],GAMMA)-K*pow(rho[lym],GAMMA))/rho[l]/(Ymed(j+1)-Ymed(j-1))*Ymed(j);
	v1[l] = sqrt(vphis);
	v1[l] -= OMEGAFRAME*r*sin(Zmed(k));
      }
    }
  }
}
