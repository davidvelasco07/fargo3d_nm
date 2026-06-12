//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void thermal_diff_b (real dt) {
//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Density);
  INPUT(DivRho); //Used to store the Laplacian
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* e   = Energy->field_cpu;
  real* rho = Density->field_cpu;
  real* lapl = DivRho->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY-2;
  int size_z = Nz+2*NGHZ-2;
  real dx = Dx;
//<\EXTERNAL>
  
//<INTERNAL>
  int i; //Reserved variables 
  int j; //for the topology
  int k; //of the kernels
  int ll;
#ifdef X
  int llxm;
  int llxp;
#endif
#ifdef Y
  int llym;
  int llyp;
#endif
  real chi;
  real lapler;
//<\INTERNAL>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=2; k<size_z; k++) {
#endif
#ifdef Y
    for(j=2; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>
#ifdef X
	ll = l;
#endif
	e[ll] += dt*lapl[ll]*rho[ll];
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
