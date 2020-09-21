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
#ifdef X
  INPUT(Mpx);
  INPUT(Vx_temp);
  OUTPUT(Vx_temp);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
#ifdef X
  real* vx = Vx_temp->field_cpu;
  real* cx = Mpx->field_cpu;
#endif
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
  real dst;
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

	
	gamma_k = alpha[n]*dt;

	if (fluidtype == GAS) {
#ifdef X
	  vx[ll] = vx[ll] + cx[ll]/( 1. + c[ll] );
#endif
	}
	else{
#ifdef X
	  //Warning, we need also the velocity of the gas here....
	  vx[ll] = s_k*cx[ll]/(1. + dt*0.5*(c[ll]+c[lxm]) ) + vx[ll]/(1.+ dt*gamma_k);
#endif
	}
	
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
