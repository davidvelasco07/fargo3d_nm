#include "fargo3d.h"

void CondInit() {

  int i,j,k;
  real vphi_t;
  real r,omega,h;
  real b0;

  real *rho   = Density->field_cpu;
  real *cs    = Energy->field_cpu;
  real *vphi  = Vx->field_cpu;
  real *vr    = Vy->field_cpu;
  real *bx    = Bx->field_cpu;
  real *by    = By->field_cpu;
  real *bz    = Bz->field_cpu;

  rho = Density->field_cpu;
  //  b0 = sqrt(2.0*SIGMA0*ASPECTRATIO*MU0/(BETA*sqrt(2.0*M_PI)));
  b0 = ASPECTRATIO*sqrt(2.0*SIGMA0*MU0/BETA);

  for (k = 0; k<Nz+2*NGHZ; k++) {
    for (j = 0; j<Ny+2*NGHY; j++) {
      for (i = 0; i<Nx; i++) {
	r = Ymed(j);
	omega = sqrt(G*MSTAR/r/r/r);
	rho[l] = SIGMA0*pow((r/R0),-SIGMASLOPE)*(1.0+NOISE*(drand48()-.5));
	cs[l] = ASPECTRATIO * pow(r/R0, FLARINGINDEX) * sqrt(G*MSTAR/r);
#ifdef ADIABATIC
	cs[l] = rho[l]*rho[l]*d[l]/(GAMMA-1.0);
#endif
	h = cs[l]/(omega*r);

	vphi_t  = omega*r*sqrt(1.0+(-SIGMASLOPE + (2.0*FLARINGINDEX-1.0) + 2.0*(MAGNETICSLOPE+1)/BETA)*h*h);
	vphi_t -= OMEGAFRAME*r;
	vphi[l] = vphi_t*(1.+ASPECTRATIO*NOISE*(drand48()-.5));
	vr[l]   = r*omega*ASPECTRATIO*NOISE*(drand48()-.5);
	//	bx[l] = cs[l]*sqrt(2.0*rho[l]*MU0/BETA);
	bx[l] = b0*pow((r/R0),MAGNETICSLOPE)*(1.0+NOISE*(drand48()-.5));
	by[l] = 0.0;
	bz[l] = 0.0;
      }
    }
  }
}
