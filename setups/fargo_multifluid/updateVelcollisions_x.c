//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void UpdateVelcollisions_x_cpu (real dt) {

//<USER_DEFINED>
  INPUT(Slope);
  INPUT(Mpx);
  INPUT(Qs);
  INPUT(Vx_temp);
  OUTPUT(Vx_temp);
//<\USER_DEFINED>

//<EXTERNAL>
  real* vx = Vx_temp->field_cpu;
  real* cx = Mpx->field_cpu;
  real* pref = Qs->field_cpu;
  real* c    = Slope->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int  n     = FluidIndex;
  int fluidtype = Fluidtype;
  real *alpha = Alpha;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real gamma_k;
  real s_k;
//<\INTERNAL>
  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=1; k<size_z; k++) {
#endif
#ifdef Y
    for (j=1; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;

	gamma_k = 0.5*(pref[ll]+pref[lxm]);
	s_k     = dt*gamma_k/(1+dt*gamma_k);

	if (fluidtype == GAS)  vx[ll] = cx[ll]/( 1. + 0.5*(c[ll]+c[lxm]) );
	else                   vx[ll] = s_k*cx[ll]/(1. + 0.5*(c[ll]+c[lxm]) ) + vx[ll]/(1.+ dt*gamma_k);
	
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
