#include "fargo3d.h"

void CondInit() {
  int i,j,k;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  real ratio;
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
  
  ratio = (PRESSURERATIO*(GAMMA+1.0)+(GAMMA-1.0))/(PRESSURERATIO*(GAMMA-1.0)+(GAMMA+1.0));
  for (i = 0; i<dim1; i++) {
    rho[i] = 1.0;
    v1[i]  = 1.0;
    e[i] = .5*rho[i]*v1[i]*v1[i]*(1./ratio-ratio)/GAMMA/(ratio-PRESSURERATIO);
    if (Q1 > 0.5) {
      rho[i] *= ratio;
      v1[i] /= ratio;
      e[i] *= PRESSURERATIO;
    }
  }
}
