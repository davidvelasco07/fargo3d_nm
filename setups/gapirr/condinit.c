#include "fargo3d.h"

void CondInit() {
  int i,j,k;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *er;
  real *rho;
  real h;
  
  real omega;
  real r, r3, temperature;

  rho = Density->field_cpu;
  e   = Energy->field_cpu;
#ifdef RADIATIVE_TRANSFER
  er  = Energyrad->field_cpu;
#endif
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      h = ASPECTRATIO*Ymed(j);
      r = Ymed(j);
      r3 = r*r*r;
      omega = sqrt(G*MSTAR/(r3));
      for (i=NGHX; i<Nx+NGHX; i++) {
	v2[l] = v3[l] = 0.0;
	v1[l] = omega*r;
#ifdef CYLINDRICAL
	rho[l] = SIGMA0*pow(r/R0,-SIGMASLOPE)*exp(-pow(Zmed(k)/h,2.0)/2.0)/(ZMAX-ZMIN);
#else
	real xi = SIGMASLOPE+1.+FLARINGINDEX;
	real beta = 1.-2*FLARINGINDEX;
	real h = ASPECTRATIO*pow(r/R0,FLARINGINDEX);
	if (FLARINGINDEX == 0.0) {
	  rho[l] = SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-beta-xi+1./(h*h));
	} else {
	  rho[l] = SIGMA0/sqrt(2.0*M_PI)/(R0*ASPECTRATIO)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-xi-beta)*					\
	    exp((1.-pow(sin(Zmed(k)),-2.*FLARINGINDEX))/2./FLARINGINDEX/(h*h));
	}
	v1[l] *= sqrt(pow(sin(Zmed(k)),-2.*FLARINGINDEX)-(beta+xi)*h*h);
	v1[l] -= OMEGAFRAME*r*sin(Zmed(k));
#endif
	e[l] = rho[l]*h*h*G*MSTAR/r/(GAMMA-1.0);
	temperature = (GAMMA-1.0)*e[l]/(rho[l]*R_MU);
#ifdef RADIATIVE_TRANSFER 
	er[l] = 4.0*STEFANK*pow(temperature,4.0)/C0; // Initial estimate of radiative energy
#endif
      }
    }
  }
}
