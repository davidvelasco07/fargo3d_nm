#include "fargo3d.h"

void InProcessRadialScan_cpu (Field *Array) {
  int i,j,k;
  real* array = Array->field_cpu;

  INPUT (Array);
  OUTPUT (Array);

#ifdef Z
  for (k = 0; k < 2*NGHZ+Nz; k++) {
#endif
#ifdef Y
    for (j = 0; j < Ny+2*NGHY; j++) {
#endif
#ifdef X
      for (i = 0; i < Nx+2*NGHX; i++) {
#endif
	if (j < NGHY)
	  array[l] = 0.0;
	else
	  array[l] += array [lym];
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
}
