#include "fargo3d.h"

//Taken from arXiv:1111.1764v2 (10 May 2012)
//This file is well defined for XY,XZ,YZ runs.

void CondInit() {
  int i,j;
  int index;
  int dim1, dim2;
  real rhom, um, L, A;

  real* rho = Density->field_cpu;
  real* e   = Energy->field_cpu;

  A = AMPLITUDE;
  L = SMOOTHINGPARAMETER;
  um = 0.5*(V1-V2);
  rhom = 0.5*(RHO1-RHO2);

#if !defined(Z)
  real* v1 = Vx->field_cpu;
  real* v2 = Vy->field_cpu;
  dim1 = Nx;
  dim2 = Ny + 2*NGHY;
#define Q1 (Xmed(j) - XMIN)/(XMAX - XMIN)
#define Q2 (Ymed(i) - YMIN)/(YMAX - YMIN)

#elif !defined(Y)
  real* v1 = Vx->field_cpu;
  real* v2 = Vz->field_cpu;
  dim1 = Nx;
  dim2 = Nz + 2*NGHZ;
#define Q1 ((Xmed(j) - XMIN)/(XMAX - XMIN))
#define Q2 ((Zmed(i) - ZMIN)/(ZMAX - ZMIN))

#elif !defined(X)
  v1 = Vy->field_cpu;
  v2 = Vz->field_cpu;
  real* dim1 = Ny + 2*NGHY;
  real* dim2 = Nz + 2*NGHZ;
#define Q1 (Ymed(j) - YMIN)/(YMAX - YMIN)
#define Q2 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)

#endif  


  for (i = 0; i<dim2; i++) {
    for (j = 0; j<dim1; j++) {
      index = j+i*dim1;
      if(Q2<0.25 && Q2>=0.0) {
	rho[index] = RHO1-rhom*exp((Q2-0.25)/L);
	v1[index] = V1-um*exp((Q2-0.25)/L);
      }
      if(Q2<0.5 && Q2>=0.25) {
	rho[index] = RHO2+rhom*exp((-Q2+0.25)/L);
	v1[index]  = V2+um*exp((-Q2+0.25)/L);
      }
      if(Q2<0.75 && Q2>=0.5) {
	rho[index] = RHO2+rhom*exp((Q2-0.75)/L);
	v1[index] = V2+um*exp((Q2-0.75)/L);
      }
      if(Q2<1.0 && Q2>=0.75) {
	rho[index] = RHO1-rhom*exp((-Q2+0.75)/L);
	v1[index] = V1-um*exp((-Q2+0.75)/L);
      }
      v2[index] = A*sin(4.0*M_PI*Q1);
      e[index] = 2.5/(GAMMA-1.);
    }
  }
}
