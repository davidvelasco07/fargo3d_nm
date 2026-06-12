#include "fargo3d.h"

void CondInit() {

  int i,j,k;

  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;

  real* vy = Vy->field_cpu;
  real* vz = Vz->field_cpu;

#ifdef MHD
  real* by = By->field_cpu;
  real* bz = Bz->field_cpu;
#endif


  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j = 0; j<Ny+2*NGHY; j++) {
      for (i = 0; i<Nx; i++) {
	e[l]   = 1.0;
	rho[l] = 1.0;
#ifdef MHD
	bz[l] = 0.0;
	by[l] = 2.0*e[l]*e[l]*rho[l]/BETA;       // Azimutal B
#endif	
	vy[l]  = 2.0*A*zmed(k);                  // Azimutal V
	vz[l]  = B*ymed(j)/(pow(ymed(j),4)+1.0); // Radial V
      }
    }
  }
}
