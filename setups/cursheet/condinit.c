#include "fargo3d.h"

void CondInit(){

  real *rho;
  real *cs;
  real *v1;
  real *v2;
  real *b1;
  real *b2;
  real Q1, Q2;
  int i,j,k;

  rho = Density->field_cpu;
  cs  = Energy->field_cpu;

  if(Nz == 1) {
    b1 = Bx->field_cpu;
    b2 = By->field_cpu; 
    v1 = Vx->field_cpu;
    v2 = Vy->field_cpu;
  }
  if(Ny == 1) {
    b1 = Bx->field_cpu;
    b2 = Bz->field_cpu; 
    v1 = Vx->field_cpu;
    v2 = Vz->field_cpu;
  }  
  if(Nx == 1) {
    b1 = By->field_cpu;
    b2 = Bz->field_cpu; 
    v1 = Vy->field_cpu;
    v2 = Vz->field_cpu;
  }
  
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx; i++) {

	if(Nz == 1) {
	  Q1 = Xmed(i);
	  Q2 = Ymed(j);
	}
	if(Ny == 1) {
	  Q1 = Xmed(i);
	  Q2 = Zmed(k);
	}
	if(Nx == 1) {
	  Q2 = Zmed(k);
	  Q1 = Ymed(j);
	}

	if(fabs(Q2)>0.5)
	  b1[l] = 1.0/sqrt(4.0*M_PI);
	else 
	  b1[l] = -1.0/sqrt(4.0*M_PI);
	b2[l] = 0.0;
	v2[l] = AMPLITUDE*sin(2.0*M_PI*Q1);
	v1[l] = 0.0;
	rho[l] = 1.0;
	cs[l]=BETA/(2.0*(GAMMA-1.0));

      }
    }
  }
}
