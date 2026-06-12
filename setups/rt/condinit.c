#include "fargo3d.h"

void CondInit() {
  
  int i,j,k;
  real* rho = Density->field_cpu;
  real* e   = Energy->field_cpu;
  real* vx  = Vx->field_cpu;
  real* vy  = Vy->field_cpu;
  real e0;
  real c;
  real A; //Amplitude
  static int init = 0;

  e0 = (2.5)/(GAMMA-1.0);

  for(k=0;k<Nz+2*NGHZ;k++) {
    for(j=0;j<Ny+2*NGHY;j++) {
      for(i=0;i<Nx;i++) {
	A = -0.01+drand48()*0.02;
	if(Ymin(j)<0.0) {
	  rho[l] = RHO1;
	  e[l] = GRAVITY*RHO1*Ymed(j)/(GAMMA-1.0) + e0;
	}
	else {
	  if (init == 0) {
	    YLAST = Ymin(j);
	    init = 1;
	  }
	  rho[l] = RHO2;
	  e[l] = GRAVITY*RHO2*(Ymed(j)-YLAST)/(GAMMA-1.0) + e0;
	}
	vx[l] = 0.0;

        vy[l] = 0*0.01/4.0 * ((1+cos(2*M_PI*Xmed(i)/0.5)) *
			    (1+cos(2*M_PI*Ymin(j)/1.5))); //Monomode
        //vy[l] = 0.01/4.0 * ((1+cos(2*M_PI*Xmed(i)/0.5))); //Monomode

//	if(Ymin(j)>-0.05 && Ymin(j)<0.05)
//	  vy[l] = 0.01/4.0 * ((1+cos(2*M_PI*Xmed(i)/0.5))); //Monomode
////	  vy[l] = A*(1+cos(8*M_PI*Xmin(i)/3.0))/2.0; //Multimode
//	else
//	vy[l] = 0.0;
      }
    }
  }
}
