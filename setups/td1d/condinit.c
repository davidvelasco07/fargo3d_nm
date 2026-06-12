#include "fargo3d.h"

void CondInit() {
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  real thdiff = THERMALDIFFUSION;
  int index;
  real* v1 = Vx->field_cpu;
  int dim1 = Nx + 2*NGHX;
  real med = .5*(XMAX + XMIN);
  for (i = 0; i<dim1; i++) {
    //e[i]   = 1.0/(GAMMA-1.0);
    rho[i] = 1.0;
    e[i] = 1.0 + 1.E-8 * exp(-.5*pow((Xmed(i)-med)/((XMAX-XMIN)/30.),2));
    v1[i]  = 0.0;
  }
}
