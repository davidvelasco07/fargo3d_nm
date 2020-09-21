//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ComputeCBcollisions_cx_cpu (real dt) {
 
//<USER_DEFINED>
  INPUT(Density);
  INPUT(Total_Density);
#ifdef X
  INPUT(Vx_temp);
  OUTPUT(Mpx);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens     = Density->field_cpu;
  real* dens_gas = Total_Density->field_cpu;
#ifdef X
  real* vx = Vx_temp->field_cpu;
  real* cx = Mpx->field_cpu;
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int n      = FluidIndex;
  int fluidtype = Fluidtype;
  real* alpha = Alpha;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real gamma_k;
  real s_k;
  real dst;
  real _cx;
  real epsilon;
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
	
	epsilon = (dens[ll]+dens[lxm])/(dens_gas[ll]+dens_gas[lxm]);
	gamma_k = alpha[n]*dt;	
	
#ifdef X	
	if (fluidtype == GAS)  _cx = 0.;
	else _cx  = epsilon*gamma_k/(1+gamma_k)*vx[l];
	
	cx[ll] += _cx;
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
