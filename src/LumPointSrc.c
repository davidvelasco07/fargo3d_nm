//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void LumPointSrc_cpu(real dt) {
  
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
  real dx = Dx;
  real luminosity = LUMOS;
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
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
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
  for(k=k0; k<size_z+k0; k++) {
#endif
#ifdef Y
    for(j=j0; j<size_y+j0; j++) {
#endif
#ifdef X
      for(i=i0; i<size_x+i0; i++) {
#endif
//<#>
	frac=1;
	if(fabs(xbody-XC)<zone_size_x(j,k) && fabs(ybody-YC)<zone_size_y(j,k) && fabs(zbody-ZC)<zone_size_z(j,k)){
#ifdef X
	  frac *= 1-fabs((xbody-XC)/zone_size_x(j,k));
#endif
#ifdef Y
	  frac *= 1-fabs((ybody-YC)/zone_size_y(j,k));
#endif
#ifdef Z
	  frac *= 1-fabs((zbody-ZC)/zone_size_z(j,k));
#endif
	  e[l] += dt*luminosity*InvVol(j, k)*frac;
	}

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

