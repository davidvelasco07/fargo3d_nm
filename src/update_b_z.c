//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Update_b_Z_cpu(Field *Q, Field *Flux) {

//<USER_DEFINED>
  INPUT(Q);
  INPUT(Flux);
  OUTPUT(Q);
//<\USER_DEFINED>

//<EXTERNAL>
  real* flux = Flux->field_cpu;
  real* qb = Q->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ-1;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Reserved variables 
  int j; //for the topology
  int k; //of the kernels
  int ll;
  int llzp;
//<\INTERNAL>

//<CONSTANT>
// real Sxj(Ny+2*NGHY);
// real Syj(Ny+2*NGHY);
// real Szj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syk(Nz+2*NGHZ);
// real Szk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
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
	ll = l;
	llzp = lzp;

	qb[ll] += ((flux[ll] - flux[llzp]) * InvVol(j,k));
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
