//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void SubStep4_b_cpu (real dt) {

//<USER_DEFINED>
#ifdef X     
  INPUT(Mmx);
#endif
#ifdef Y
  INPUT(Mmy);
#endif
#ifdef Z
  INPUT(Mmz);
#endif
  INPUT(Density);
  INPUT(Energy);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
#ifdef X
  real* fx = Mmx->field_cpu;
#endif
#ifdef Y
  real* fy = Mmy->field_cpu;
#endif
#ifdef Z
  real* fz = Mmz->field_cpu;
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX-1; 
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
  real constant = dt*R_MU/(GAMMA-1.0);
  real rho0 = RHO0;
//<\EXTERNAL>
  
//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
  real div;
  real r;
  real h;
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
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=1; j<size_y; j++) {
#endif
#ifdef X
      for(i=1; i<size_x; i++) {
#endif
//<#>
	ll = l;
	div = 0.0;

#ifdef X
	div += (fx[lxp]-fx[ll]);
#endif
#ifdef Y
	div += (fy[lyp]-fy[ll]);
#endif
#ifdef Z
	div += (fz[lzp]-fz[ll]);
#endif
	
#if THDIFFUSION>2//Planet-disc cases	
	e[ll] -= constant * div *  InvVol(j,k);
#else	
	e[ll] -= constant * rho0 * div *  InvVol(j,k);
#endif
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
