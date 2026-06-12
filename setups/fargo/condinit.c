#include "fargo3d.h"

void CondInit() {
  
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Pot);
  OUTPUT(Vx);
  OUTPUT(Vy);

  int i,j,k;
  real r, omega;
  
  real *vphi = Vx->field_cpu;
  real *vr   = Vy->field_cpu;
  real *rho  = Density->field_cpu;
#ifdef ADIABATIC
  real *e   = Energy->field_cpu;
#endif
#ifdef ISOTHERMAL
  real *cs   = Energy->field_cpu;
#endif

  i = j = k = 0;
  
  for (j=0; j<Ny+2*NGHY; j++) {
    for (i=0; i<Nx+2*NGHX; i++) {
      
      r = Ymed(j);
      omega = sqrt(G*MSTAR/r/r/r);
      
      rho[l] = SIGMA0*pow(r/R0,-SIGMASLOPE)*(1.0+NOISE*(drand48()-.5));
      cs[l]  = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*omega*r;
      
#ifdef ADIABATIC
      e[l] = cs[l]*cs[l]*rho[l]/(GAMMA-1.0);
#endif
      
      vphi[l] = omega*r*sqrt(1.0+pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
			(2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));
      vphi[l] -= OMEGAFRAME*r;
      
      vphi[l] *= (1.+ASPECTRATIO*NOISE*(drand48()-.5));
      vr[l] = cs[l]*NOISE*(drand48()-.5);
    }
  } 
}
