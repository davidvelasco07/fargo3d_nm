#include "fargo3d.h"
#define EPS 0.0001

void CondInit() {
  int i,j,k;
  real kernel;
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  int index;
  real projx, projy;
  
#ifndef Z
  real* v1 = Vx->field_cpu;
  real* v2 = Vy->field_cpu;
  int dim1 = Nx + 2*NGHX;
  int dim2 = Ny + 2*NGHY;
  projx = 1./sqrt(1.+pow((YMAX-YMIN)/(XMAX-XMIN),2.0));
  projy = 1./sqrt(1.+pow((XMAX-XMIN)/(YMAX-YMIN),2.0));
#define Q1 (Xmed(j) - XMIN)/(XMAX - XMIN)
#define Q2 (Ymed(i) - YMIN)/(YMAX - YMIN)
#endif
#ifndef Y
  real* v1 = Vx->field_cpu;
  real* v2 = Vz->field_cpu;
  int dim1 = Nx + 2*NGHX;
  int dim2 = Nz + 2*NGHZ;
  projx = 1./sqrt(1.+pow((ZMAX-ZMIN)/(XMAX-XMIN),2.0));
  projy = 1./sqrt(1.+pow((XMAX-XMIN)/(ZMAX-ZMIN),2.0));
#define Q1 (Xmed(j) - XMIN)/(XMAX - XMIN)
#define Q2 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#endif
#ifndef X
  real* v1 = Vy->field_cpu;
  real* v2 = Vz->field_cpu;
  int dim1 = Ny + 2*NGHY;
  int dim2 = Nz + 2*NGHZ;
  projx = 1./sqrt(1.+pow((ZMAX-ZMIN)/(YMAX-YMIN),2.0));
  projy = 1./sqrt(1.+pow((YMAX-YMIN)/(ZMAX-ZMIN),2.0));
#define Q1 (Ymed(j) - YMIN)/(YMAX - YMIN)
#define Q2 (Zmed(i) - ZMIN)/(ZMAX - ZMIN)
#endif
  
  for (i = 0; i<dim2; i++) {
    for (j = 0; j<dim1; j++) {
      index = j+i*dim1;
      
      kernel = exp(-(((Q1-.5)*(Q1-.5)+(Q2-.5)*(Q2-.5))*25.));
      rho[index] = (1.+kernel*EPS*cos(20.*M_PI*(Q1+Q2)));
      v1[index] = CS*EPS*cos(20.*M_PI*(Q1+Q2))*kernel*projx;
      v2[index] = CS*EPS*cos(20.*M_PI*(Q1+Q2))*kernel*projy;
      e[index] = CS;
#ifdef ADIABATIC
      e[index] = (1.+GAMMA*kernel*EPS*cos(20.*M_PI*(Q1+Q2)))*CS*CS/(GAMMA*(GAMMA-1.0));
#endif
      
    }
  }
}
