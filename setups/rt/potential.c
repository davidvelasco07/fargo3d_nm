//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void compute_potential(real dt) {
  static init = 0;
  
  if (init == 0) {
    FARGO_SAFE(Potential()); // Gravitational potential from star and planet(s)
    init = 1;
  }
  return;
}

void Potential_cpu() {
  
//<USER_DEFINED>
  OUTPUT(Pot);
//<\USER_DEFINED>
  
//<EXTERNAL>
  real* pot  = Pot->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
//<\INTERNAL>

//<CONSTANT>
// real ymin(Ny+2*NGHY+1);
// real GRAVITY(1);
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
      for (i=0; i<size_x; i++) {
#endif
//<#>
	pot[l] =  -GRAVITY*ymed(j);
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
