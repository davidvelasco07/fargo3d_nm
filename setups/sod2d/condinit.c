#include "fargo3d.h"

void CondInit() {
  int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  int index;
#ifndef Z
  real* v1 = Vx->field_cpu;
  real* v2 = Vy->field_cpu;
  int dim1 = Nx;
  int dim2 = Ny + 2*NGHY;
#define Q1 (Xmed(j) - XMIN)/(XMAX - XMIN)
#define Q2 (Ymed(i) - YMIN)/(YMAX - YMIN)
#endif
#ifndef Y
  real* v1 = Vx->field_cpu;
  real* v2 = Vz->field_cpu;
  int dim1 = Nx;
  int dim2 = Nz + 2*NGHZ;
#define Q1 (Xmed(j) - XMIN)/(XMAX - XMIN)
#define Q2 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#endif
#ifndef X
  real* v1 = Vy->field_cpu;
  real* v2 = Vz->field_cpu;
  int dim1 = Ny + 2*NGHY;
  int dim2 = Nz + 2*NGHZ;
#define Q1 (Ymed(j) - YMIN)/(YMAX - YMIN)
#define Q2 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#endif

  for (i = 0; i<dim2; i++) {
    for (j = 0; j<dim1; j++) {
      index = j+i*dim1;
      
      e[index]   = 1.0/(GAMMA-1.0); //Overwritten in what follows if adiabatic
      rho[index] = 1.0;
      v1[index]  = 0.0;
      v2[index]  = 0.0;
      
      if ((Q1+Q2/.6) > 1.) {
      	rho[index] = 0.125;
#ifdef ADIABATIC
	e[index] = 0.1/(GAMMA-1.0);
#endif
      }
    }
  }
}
