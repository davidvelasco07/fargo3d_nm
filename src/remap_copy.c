//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ReMapCopy_cpu() {

  //Temporal function, should be replaced with a propper copy function

//<USER_DEFINED>
  INPUT(Pressure);
  INPUT(Mpx);
  INPUT(Mmx);
  INPUT(Mmy);
  OUTPUT(Density);
  OUTPUT(Vx);
  OUTPUT(Vy);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  int size_x = Nx+2*NGHX;
  int size_y= Ny+2*NGHY;
  int size_z= Nz+2*NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real* rho = Density->field_cpu;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
  real* cs = Energy->field_cpu;
  real* rho_temp = Pressure->field_cpu;
  real* vx_temp = Mpx->field_cpu;
  real* vy_temp = Mmx->field_cpu;
  real* cs_temp = Mmy->field_cpu;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
//<\INTERNAL>

  i = j = k = 0;

//<MAIN_LOOP>
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
	rho[ll] = rho_temp[ll];
	vx[ll] = vx_temp[ll];
	vy[ll] = vy_temp[ll];
	cs[ll] = cs_temp[ll];

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
