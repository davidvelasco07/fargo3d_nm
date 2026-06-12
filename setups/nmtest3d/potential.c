//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void compute_potential(real dt) {
  FARGO_SAFE(Potential());
}

void Potential_cpu() {
  
//<USER_DEFINED>
  OUTPUT(Pot);
  real planetmass_taper;
  if (MASSTAPER == 0.0)
    planetmass_taper = 1.0;
  else
    planetmass_taper = (PhysicalTime >= MASSTAPER ? 1.0 : .5*(1.0-cos(M_PI*PhysicalTime/MASSTAPER)));
//<\USER_DEFINED>

//<EXTERNAL>
  real* pot  = Pot->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real xbody = XBODY;
  real ybody = YBODY;
  real zbody = ZBODY;
  real mass  = MASSBODY*planetmass_taper;
  real sm2   = SMOOTHING*SMOOTHING;
//<\EXTERNAL>
  
//<INTERNAL>
  int i;
  int j;
  int k;
  int n;
  real dist2;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
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
	dist2 = (XC-xbody)*(XC-xbody)+		\
		(YC-ybody)*(YC-ybody)+		\
		(ZC-zbody)*(ZC-zbody);
	pot[l] = -G*mass/sqrt(dist2+sm2);
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
//<\MAIN_LOOP>
#endif
}
