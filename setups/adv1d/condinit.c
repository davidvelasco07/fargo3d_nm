#include "fargo3d.h"

void CondInit() {
  int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  int index;
#ifdef Z
  real* v1 = Vz->field_cpu;
  int dim1 = Nz + 2*NGHZ;
#define Q1 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#endif
#ifdef Y
  real* v1 = Vy->field_cpu;
  int dim1 = Ny + 2*NGHY;
#define Q1 (Ymed(i) - YMIN)/(YMAX - YMIN)
#endif
#ifdef X
  real* v1 = Vx->field_cpu;
  int dim1 = Nx;
#define Q1 (Xmed(i) - XMIN)/(XMAX - XMIN)
#endif
  
  for (i = 0; i<dim1; i++) {
    e[i]   = CS;
    rho[i] = RHO0;
    v1[i]  = V0;
    if ((Q1 >= 0.25) && (Q1 < 0.75))
      rho[i] = 1.0;
  }
}
