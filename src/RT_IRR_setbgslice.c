//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void SetBackgroundSlice_cpu (Field* array) {


//<USER_DEFINED>
  INPUT2D(LocalDepth);
  INPUT(array); // <== IMPORTANT. Only one slice of 'array' is
		// modified, therefore others need to be up to date
  OUTPUT(array);
//<\USER_DEFINED>

//<EXTERNAL>
  real* f   = array->field_cpu;
  real* s   = LocalDepth->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
#ifdef __GPU
  int pitch2d = LocalDepth->pitch;
#else
  int pitch2d = Nx+2*NGHX;
#endif
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>
  
//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
//<\INTERNAL>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=0; k<size_z; k++) {
#endif
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>

	if (j == size_y-NGHY-1) // We store the last slice in radius in a 2D, ZX array.
	  f[l] = s[i+pitch2d*k];

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
