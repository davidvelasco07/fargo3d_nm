//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void LumSphericalSrc_cpu(real dt,int LumCells) {
  
//<USER_DEFINED>
  INPUT(Energy);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* e   = Energy->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real xbody = XBODY;
  real ybody = YBODY;
  real zbody = ZBODY;
  real luminosity = LUMOS;
  real sm2=SMOOTHING*SMOOTHING;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int i0 = 0;
  int j0 = 0;
  int k0 = 0;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  real frac;
  real x;
  real y;
  real z;
  real r2;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real Syk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
//<\CONSTANT>

//<MAIN_LOOP>
  
  i = j = k = 0;
#ifdef Z
  for(k=k0; k<size_z+k0; k++) {
#endif
#ifdef Y
    for(j=j0; j<size_y+j0; j++) {
#endif
#ifdef X
      for(i=i0; i<size_x+i0; i++) {
#endif
//<#>
	r2=0;
#ifdef X
	x = XC-xbody;
	r2 += x*x;
#endif
#ifdef Y
	y = YC-ybody;
	r2 += y*y;
#endif
#ifdef Z
	z = ZC-zbody;
	r2 += z*z;
#endif	
	frac = (r2<=sm2) ? dt*luminosity/LumCells:0;
	e[l] += InvVol(j, k)*frac;
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

