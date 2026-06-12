#include "fargo3d.h"

void CondInit() {
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);
  int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  int index;
  real* v1 = Vx->field_cpu;
  real* v2 = Vy->field_cpu;
  int dim1 = Nx + 2*NGHX;
  int dim2 = Ny + 2*NGHY;
  real medx = .5*(YMAX + YMIN);
  real medy = 0.0;
  real x, y;
  for (j = 0; j<dim2; j++) {
    for (i = 0; i<dim1; i++) {
      index = i+j*dim1;
      rho[index] = 1.0;
      x = Ymed(j)*cos(Xmed(i));
      y = Ymed(j)*sin(Xmed(i));
      e[index]   = 1.0 + 1.E-8 * exp(-.5 * (pow((x-medx)/(0.1),2) + pow((y-medy)/(0.1),2))); 
      v1[index]  = 0.0;
      v2[index]  = 0.0;
      
    }
  }
}
