//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RT_IRR_TauStellar_cpu () {
//<USER_DEFINED>
  INPUT(Density);
  OUTPUT(Tau);
//<\USER_DEFINED>

  /* Under the standard practice of FARGO3D, KAPPASTAR should have a
     scaling rule and should be expressed in the current unit system
     in a given parameter file. However, this would be hardly legible
     as we use dimensionless variables for most of the parameters,
     which are therefore O(1). We prefer here to enforce its
     definition in cgs, and we convert in this function on the
     fly. Harmless warnings are issued at build time. */


//<EXTERNAL>
  real* rho     = Density->field_cpu;
  real* tau     = Tau->field_cpu;
  real  opacity = KAPPASTAR*(MSTAR_CGS/MSTAR)*pow(R0/R0_CGS,2.0);
  int   pitch   = Pitch_cpu;
  int   stride  = Stride_cpu;
  int   size_x  = Nx+2*NGHX;
  int   size_y  = Ny+2*NGHY;
  int   size_z  = Nz+2*NGHZ;
//<\EXTERNAL>
  
//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
//<\INTERNAL>

//<CONSTANT>
// real ymin(Ny+2*NGHY+1);
//<\CONSTANT>

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
	ll = l;
	tau[ll] = rho[ll]*opacity*(ymin(j+1)-ymin(j));
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

//<LAST_BLOCK>
  RadialScan (Tau);
//<\LAST_BLOCK>
}
