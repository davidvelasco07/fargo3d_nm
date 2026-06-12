//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void MultiplyByRho_cpu(Field *Q, Field *Mult) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Q);
  OUTPUT(Mult);
//<\USER_DEFINED>


//<EXTERNAL>
  real* q      = Q->field_cpu;
  real* multrho = Mult->field_cpu;
  real* rho    = Density->field_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int stride = Stride_cpu;
  int pitch  = Pitch_cpu;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
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
      for (i=0; i<size_x; i++) {
#endif
//<#>
	ll = l;
	multrho[ll] = q[ll]*rho[ll];
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
