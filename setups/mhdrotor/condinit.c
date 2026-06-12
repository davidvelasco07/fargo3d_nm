#include "fargo3d.h"

real taper(real r) {
  return (R2-r)/(R2-R1);
}

void CondInit() {
  int i,j,k;
  real Q1,Q2;

  real *v1;
  real *v2;
  real *b1;
  real *b2;
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

  if(Nz == 1) {
    b1 = Bx->field_cpu;
    b2 = By->field_cpu; 
    v1 = Vx->field_cpu;
    v2 = Vy->field_cpu;
    printf("MHD ROTOR XY SETUP\n");
  }
  if(Ny == 1) {
    b1 = Bx->field_cpu;
    b2 = Bz->field_cpu; 
    v1 = Vx->field_cpu;
    v2 = Vz->field_cpu;
    printf("MHD ROTOR XZ SETUP\n");
  }
  if(Nx == 1) {
    b1 = By->field_cpu;
    b2 = Bz->field_cpu; 
    v1 = Vy->field_cpu;
    v2 = Vz->field_cpu;
    printf("MHD ROTOR YZ SETUP\n");
  }

  rho = Density->field_cpu;
  cs  = Energy->field_cpu;
  
  i = j = k = 0;
  
  real rmed, rymin, rzmin;

#ifdef Z
  for (k=0; k<Nz+2*NGHZ; k++) {
#endif
#ifdef Y
    for (j=0; j<Ny+2*NGHY; j++) {
#endif
#ifdef X
      for (i=0; i<Nx; i++) {
#endif
	if(Nz == 1) {
	  Q1 = (Xmed(i)-XMIN)/(XMAX-XMIN);
	  Q2 = (Ymed(j)-YMIN)/(YMAX-YMIN);
	}
	if(Ny == 1) {
	  Q1 = (Xmed(i)-XMIN)/(XMAX-XMIN);
	  Q2 = (Zmed(k)-ZMIN)/(ZMAX-ZMIN);
	}
	if(Nx == 1) {
	  Q1 = (Ymed(j)-YMIN)/(YMAX-YMIN);
	  Q2 = (Zmed(k)-ZMIN)/(ZMAX-ZMIN);
	}

	b1[l] = 5.0/sqrt(4.0*M_PI);
	b2[l] = 0.0;
	cs[l] = P/(GAMMA-1.0);

	rmed = sqrt((ymed(j)-0.5)*(ymed(j)-0.5)+(zmed(k)-0.5)*(zmed(k)-0.5));
	rymin = sqrt((ymin(j)-0.5)*(ymin(j)-0.5)+(zmed(k)-0.5)*(zmed(k)-0.5));
	rzmin = sqrt((ymed(j)-0.5)*(ymed(j)-0.5)+(zmin(k)-0.5)*(zmin(k)-0.5));

	rho[l] = 1.0+9.0*taper(rmed);
	v1[l] = -taper(rymin)*V0*(zmed(k)-0.5)/R1;
	v2[l] = taper(rzmin)*V0*(ymed(j)-0.5)/R1;
	
	if (rmed<R1)
	  rho[l] = 10.0;
	if (rymin<R1)
	  v1[l] = -V0*(zmed(k)-0.5)/R1;
	if (rzmin<R1)
	  v2[l] = V0*(ymed(j)-0.5)/R1;

	if (rmed>R2 || rmed<0.05)
	  rho[l] = 1.0;
	if (rymin>R2 || rymin<0.05)
	  v1[l] = 0.0;
	if (rzmin>R2 || rzmin<0.05)
	  v2[l] = 0.0;
	
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
