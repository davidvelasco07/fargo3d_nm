//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void OffsetPatches_cpu (Field *F) {

//<USER_DEFINED>
  INPUT(F);
  INPUT2D(LocalDepth);
  OUTPUT(F);
//<\USER_DEFINED>

//<EXTERNAL>
  real* f    = F->field_cpu;
  real* s    = LocalDepth->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
#ifdef __GPU
  int pitch2d = LocalDepth->pitch;
#else
  int pitch2d = Nx+2*NGHX;
#endif
  int ny     = Ny;
//<\EXTERNAL>
  
//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int l2d; //NOT TO BE CONFUSED WITH l2D WHICH IS A PREPROCESSOR DEFINITION
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
	l2d = i + pitch2d*k;
	if (j != ny+NGHY-1) //We want to avoid race conditions
	  f[l] += s[l2d] - f[i+(ny+NGHY-1)*pitch+k*stride]; //CPU&GPU version
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
