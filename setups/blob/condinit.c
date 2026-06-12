#include "fargo3d.h"

void CondInit() {

  int i,j,k;

  real* rho = Density->field_cpu;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
  real* e  = Energy->field_cpu;

  int dim1 = Nx;
  int dim2 = Ny + 2*NGHY;

  int index;

#define Q1 (Xmed(j) - XMIN)/(YMAX - YMIN)
#define Q2 (Ymed(i) - YMIN)/(YMAX - YMIN)

  for (i = 0; i<dim2; i++) {
    for (j = 0; j<dim1; j++) {
      index = j+i*dim1;

      rho[index] = RHO1;
      vx[index] = V1;
      vy[index] = 0.0;
      e[index]  = CS*CS/(GAMMA-1.0);
      
      if(sqrt((Q1-XBLOB)*(Q1-XBLOB)+(Q2-YBLOB)*(Q2-YBLOB)) < RBLOB) {
	rho[index] = RHO2;
	vx[index] = V2;
      }
    }
  }
}
