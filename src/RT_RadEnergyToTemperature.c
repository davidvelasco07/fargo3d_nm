//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void RadEnergyToTemperature_cpu () {

//<USER_DEFINED>
  INPUT(Eta1);
  INPUT(Eta2);
  INPUT(Energyrad);
  OUTPUT(Temperature);
//<\USER_DEFINED>


//<EXTERNAL>
  real* energyrad = Energyrad->field_cpu;
  real* temper    = Temperature->field_cpu;
  real* eta1      = Eta1->field_cpu;
  real* eta2      = Eta2->field_cpu;
  int pitch       = Pitch_cpu;
  int stride      = Stride_cpu;
  int size_x      = Nx+NGHX;
  int size_y      = Ny+NGHY;
  int size_z      = Nz+NGHZ;
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
	temper[ll] = eta1[ll] + energyrad[ll]*eta2[ll];
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
