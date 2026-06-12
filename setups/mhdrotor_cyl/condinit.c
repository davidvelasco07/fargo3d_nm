#include "fargo3d.h"

real taper(real r) {
  return (R2-r)/(R2-R1);
}

real pot(real r, real phi) {
  //Potential vector A=-r*cos(phi)*B0. --> B = B0 * \vec{e}_j
  real B0, A;
  B0 = 5.0/sqrt(4.0*M_PI);
  A = -r*cos(phi)*B0;	 
  return A;
}

void CondInit() {
  int i,j,k;
  real Q1,Q2;

  real *v1;
  real *v2;
  real *v3;
  real *b1;
  real *b2;
  real *b3;
  real *rho;
  real *cs;

  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);
  OUTPUT(Vz);
  OUTPUT(Bx);
  OUTPUT(By);
  OUTPUT(Bz);

  b1 = Bx->field_cpu;
  b2 = By->field_cpu; 
  b3 = Bz->field_cpu; 
  v1 = Vx->field_cpu;
  v2 = Vy->field_cpu;
  v3 = Vz->field_cpu;
  printf("MHD ROTOR R-PHI SETUP\n");

  rho = Density->field_cpu;
  cs  = Energy->field_cpu;
  
  i = j = k = 0;
  
  real rmed;

#ifdef Z
  for (k=0; k<Nz+2*NGHZ-1; k++) {
#endif
#ifdef Y
    for (j=0; j<Ny+2*NGHY-1; j++) {
#endif
#ifdef X
      for (i=0; i<Nx; i++) {
#endif

	b1[l] = - (pot(Ymin(j+1),Xmin(i)) - pot(Ymin(j),Xmin(i))) / (Ymin(j+1)-Ymin(j));   //-dA/dr
	b2[l] = 1/Ymin(j)*(pot(Ymin(j),Xmin(i+1)) - pot(Ymin(j),Xmin(i)))/(Xmin(i+1)-Xmin(i)); //1/r dA/dphi
	b3[l] = 0.0;
	
	rmed = Ymed(j);
	
	rho[l] = 1.0+9.0*taper(rmed);
	v1[l] = V0*taper(rmed)*rmed/R1;
	v2[l] = 0.0;
	v3[l] = 0.0;

	if (rmed<R1){
	  rho[l] = 10.0;
	  v1[l] = V0*rmed/R1;
	}
	if (rmed>R2){
	  rho[l] = 1.0;
	  v1[l] = 0.0;
	}
	if (rmed < 0.05) {
	  rho[l] = 1.0;
	  v1[l] = 0.0;
	}

	cs[l] = P/(GAMMA-1.0);//sqrt(P/rho[l]);
	
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
}
