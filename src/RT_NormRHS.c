//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void NormRHS_cpu () {

//<USER_DEFINED>
  INPUT (RHS);
  OUTPUT (Resid);
//<\USER_DEFINED>


//<EXTERNAL>
  real* rhs    = RHS->field_cpu;
  real* resid  = Resid->field_cpu;
  int pitch    = Pitch_cpu;
  int stride   = Stride_cpu;
  int size_x   = Nx+NGHX;
  int size_y   = Ny+NGHY;
  int size_z   = Nz+NGHZ;
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
  for (k=NGHZ; k<size_z; k++) {
#endif
#ifdef Y
    for (j=NGHY; j<size_y; j++) {
#endif
#ifdef X
      for (i=NGHX; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
	resid[ll] = fabs(rhs[ll]);
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
