#include "fargo3d.h"

#define Q1 (Xmed(i) - XMIN)/(XMAX - XMIN)
#define Q2 (Ymed(j) - YMIN)/(YMAX - YMIN)
#define Q3 (Zmed(k) - ZMIN)/(ZMAX - ZMIN)

void CondInit() {

  int i,j,k;
  int dim1;
  int dim2;
  int dim3;

  real *vx;
  real *vy;
  real *vz;
  real *e;
  real *rho;

  real vol = 4./3. * M_PI * pow(RADIUS,3);

  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  vx  = Vx->field_cpu;
  vy  = Vy->field_cpu;
  vz  = Vz->field_cpu;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx; i++) {
	vx[l] = vy[l] = vz[l] = 0.0;
	rho[l] = 1.0;
	e[l] = 1.0e-6;
	if (sqrt(pow((Q1-0.5),2) + 
		 pow((Q2-0.5),2) + 
		 pow((Q3-0.5),2)) < RADIUS)
	  e[l] *= E0/vol;
      }
    }
  }
}
