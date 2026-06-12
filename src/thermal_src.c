//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ThermalSrc_cpu(real dt) {
  
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
#ifdef CARTESIAN
  int size_x = 2;
  int isrc = floor((XBODY-Xmin(0))/(Xmin(1)-Xmin(0))-0.5);
  int size_y = 2;
  int jsrc = floor((YBODY-Ymin(0))/(Ymin(1)-Ymin(0))-0.5);
  int size_z = 2;
  int ksrc = floor((ZBODY-Zmin(0))/(Zmin(1)-Zmin(0))-0.5);
#elif CYLINDRICAL
  int size_x = Nx;
  int isrc = NGHX;
  int size_y = 1;
  int jsrc = NGHY;
  int size_z = 1;
  int ksrc =  floor((ZBODY-Zmin(0))/(Zmin(1)-Zmin(0))-0.5);
#else
  int size_x = 2;
  int isrc = floor((XBODY-Xmin(0))/(Xmin(1)-Xmin(0))-0.5);
  int size_y = 2;
  int jsrc = floor((YBODY-Ymin(0))/(Ymin(1)-Ymin(0))-0.5);
  int size_z = 2;
  int ksrc = floor((ZBODY-Zmin(0))/(Zmin(1)-Zmin(0))-0.5);
#endif
  real luminosity = LUMOS;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  real fracx;
  real fracy;
  real fracz;
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

#ifdef CARTESIAN
  fracz = fracy = fracx = 0;
#endif
#ifdef CYLINDRICAL
  fracz = fracy = 1;
  fracx = 1./Nx;
#endif

#ifdef Z
  for(k=ksrc; k<size_z+ksrc; k++) {
#endif
#ifdef Y
    for(j=jsrc; j<size_y+jsrc; j++) {
#endif
#ifdef X
      for(i=isrc; i<size_x+isrc; i++) {
#endif
//<#>
#ifdef CARTESIAN
#ifdef X
	fracx = 1-fabs((xbody-xmed(i))/(xmed(i+1)-xmed(i)));
#endif
#ifdef Y
	fracy = 1-fabs((ybody-ymed(j))/(ymed(j+1)-ymed(j)));
#endif
#ifdef Z
	fracz = 1-fabs((zbody-zmed(k))/(zmed(k+1)-zmed(k)));
#endif	
#endif
	e[l] += dt*luminosity*InvVol(j, k)*(fracx*fracy*fracz);
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

