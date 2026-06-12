#include "fargo3d.h"

void CondInit() {

  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);
  OUTPUT(Vz);

  int i,j,k;

  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  real* v1 = Vx->field_cpu;
  real* v2 = Vy->field_cpu;
  real* v3 = Vz->field_cpu;
  
  int dim1 = Nx;
  int dim2 = Ny+2*NGHY;
  int dim3 = Nz+2*NGHZ;

#define Q1 (Xmed(i) - XMIN)/(XMAX - XMIN)
#define Q2 (Ymed(j) - YMIN)/(YMAX - YMIN)
#define Q3 (Zmed(k) - ZMIN)/(ZMAX - ZMIN)

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx; i++) {
	v1[l] = v2[l] = v3[l] = 0.0;
	rho[l] = 10.0/(YMAX-YMIN);
#ifdef ADIABATIC
	e[l] = 1.0/(GAMMA-1.0);
#else
	e[l]=CS;
#endif
	if(Q1+Q2+Q3>1.0) {
	  rho[l] = 0.125;	  
#ifdef ADIABATIC
	  e[l] = 0.1/(GAMMA-1.0);
#endif
	}
      }
    }
  }
}
