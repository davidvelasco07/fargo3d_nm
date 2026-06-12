//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ComputeEtas_cpu (real dt) {

//<USER_DEFINED>
  INPUT(Temperature);
  INPUT(OpaPlanck);
#ifdef IRRADIATION
  INPUT(StellarRad);
#endif
  INPUT(QPlus);
  INPUT(Density);
  OUTPUT(Eta1);
  OUTPUT(Eta2);
  OUTPUT(Energyrad);
//<\USER_DEFINED>


//<EXTERNAL>
  real* temper     = Temperature->field_cpu;
  real* opap       = OpaPlanck->field_cpu;
#ifdef IRRADIATION
  real* stellarrad = StellarRad->field_cpu;
#endif
  real* qplus      = QPlus->field_cpu;
  real* dens       = Density->field_cpu;
  real* eta1       = Eta1->field_cpu;
  real* eta2       = Eta2->field_cpu;
  real* energyrad  = Energyrad->field_cpu;
  int pitch        = Pitch_cpu;
  int stride       = Stride_cpu;
  int size_x       = Nx+2*NGHX;
  int size_y       = Ny+2*NGHY;
  int size_z       = Nz+2*NGHZ;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real denom;
  real cv;
//<\INTERNAL>
  
//<CONSTANT>
// real GAMMA(1);
// real ZMIN(1);
// real ZMAX(1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
	cv = R_MU/(GAMMA-1.0);
	/* In the expression below, qplus[ll] is the gas heating
	   rate. It is the sum of the viscous heating and stellar
	   irradiation, if any. It is required to embed both kinds of
	   heating in the same variable in order to communicate only
	   one field (qplus).  */
	denom = 1.0+16.0*dt*opap[ll]*STEFANK*pow(temper[ll],3.0)/cv;
	eta1[ll] = temper[ll]+12.0*dt*opap[ll]/cv*STEFANK*pow(temper[ll],4.0)+ \
	  dt*(qplus[ll])/(dens[ll]*cv);
	eta1[ll] = eta1[ll]/denom;
	eta2[ll] = dt*C0*opap[ll]/(cv*denom);	
	energyrad[ll] = (temper[ll]-eta1[ll])/eta2[ll];
	if (zmed(k) < ZMIN) {
	  energyrad[ll] = 4.0*STEFANK*pow(TCMB,4.0)/C0;
	}
#ifndef HALFDISK
	if (zmed(k) > ZMAX) {
	  energyrad[ll] = 4.0*STEFANK*pow(TCMB,4.0)/C0;
	}
#endif

//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
