//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void DragForce_SumC_cpu (real dt) {
 
//<USER_DEFINED>
  INPUT(Density);
  INPUT(Total_Density);
  INPUT(Qs);
  OUTPUT(Slope);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens     = Density->field_cpu;
  real* dens_gas = Total_Density->field_cpu;
  real* c        = Slope->field_cpu;
  real* pref     = Qs->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real alphak;
  real sk;
  real _c;
  real epsilon;
//<\INTERNAL>

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

	epsilon = dens[ll]/dens_gas[ll];
	alphak  = pref[ll];	
	sk      = dt*alphak/(1+dt*alphak);

	if (fluidtype == GAS)  _c = 0.;
	else _c  = epsilon*sk;

	c[ll] += _c;
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
