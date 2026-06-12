#include "fargo3d.h"

real pot(real r) {
  if (r < R1)
    return A0*(R1-r);
  else
    return 0.0;
}

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
  
  for (k=0; k<Nz+2*NGHZ-1; k++) {
    for (j=0; j<Ny+2*NGHY-1; j++) {
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
	  Q2 = Zmin(k);
	  Q1 = Ymin(j);
	}

	rho[l] = 1.0;
	cs[l]=1.0/(GAMMA-1.0); //Pressure equal to 1
	v1[l] = V0X;
	v2[l] = V0Y;
	
	b1[l] = (pot(sqrt(ymin(j)*ymin(j)+zmin(k+1)*zmin(k+1))) - pot(sqrt(ymin(j)*ymin(j)+zmin(k)*zmin(k))))/(zmin(k+1)-zmin(k)); //(A0*(R1-sqrt(ymin(j)*ymin(j)+zmin(k+1)*zmin(k+1))) - 
		// A0*(R1-sqrt(ymin(j)*ymin(j)+zmin(k)*zmin(k))))/(zmin(k+1)-zmin(k));

	b2[l] = -(pot(sqrt(ymin(j+1)*ymin(j+1)+zmin(k)*zmin(k))) - pot(sqrt(ymin(j)*ymin(j)+zmin(k)*zmin(k))))/(ymin(j+1)-ymin(j));

	//	b2[l] = -(A0*(R1-sqrt(ymin(j+1)*ymin(j+1)+zmin(k)*zmin(k))) - 
	//		  A0*(R1-sqrt(ymin(j)*ymin(j)+zmin(k)*zmin(k))))/(ymin(j+1)-ymin(j));
      }
    }
  }
}
