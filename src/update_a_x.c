//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Update_a_X_cpu(real dt, Field *Qs, Field *Vx_t, Field *Flux) {

//<USER_DEFINED>
  INPUT(Qs);
  INPUT(Vx_t);
  INPUT(DensStar);
  OUTPUT(Flux);
//<\USER_DEFINED>

//<EXTERNAL>
  real* flux = Flux->field_cpu;
  real* qs = Qs->field_cpu;
  real* vx = Vx_t->field_cpu;
  real* rho_s = DensStar->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
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

	flux[ll] = vx[ll]*qs[ll]*rho_s[ll]*SurfX(j,k)*dt;
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
