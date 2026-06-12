#include "fargo3d.h"

void CondInit() {

  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Pot);
  OUTPUT(Vx);
  OUTPUT(Vy);
  OUTPUT(Vz);
#ifdef LABELED
  OUTPUT(Label);
#endif

  int i,j,k;

  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  real* v1 = Vx->field_cpu;
  real* v2 = Vy->field_cpu;
  real* v3 = Vz->field_cpu;
#ifdef LABELED
  real* label = Label->field_cpu;
#endif
  int dim1 = Nx+2*NGHX;
  int dim2 = Ny+2*NGHY;
  int dim3 = Nz+2*NGHZ;

  for (k=0; k<dim3; k++) {
    for (j=0; j<dim2; j++) {
      for (i=0; i<dim1; i++) {
	v1[l] = v2[l] = 0.0;
	v3[l] = WINDVEL;
	rho[l] = RHO0;
	e[l] = CS*CS*RHO0/(GAMMA*(GAMMA-1.0));
#ifdef LABELED
	label[l] = Zmed(k);
#endif
      }
    }
  }
}
