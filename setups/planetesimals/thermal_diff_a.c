//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void thermal_diff_a () {
//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Density);
  OUTPUT(Energy);
  OUTPUT(DivRho); //Used to store the Laplacian
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho = Density->field_cpu;
  real* e   = Energy->field_cpu;
  real* lapl = DivRho->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
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
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=1; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>
#ifdef X
	ll = l;
	llxm = lxm;
	llxp = lxp;
#endif
#ifdef Y
	llyp = lyp;
	llym = lym;
#endif
	chi = 8.0*M_PI*(GAMMA-1.0)*STEFANK*pow(e[ll]*(GAMMA-1.0)/(rho[ll]*R_MU),4.0)/ \
	  (3.0*OPACITY*rho[ll]*rho[ll])*ymed(j)*ymed(j)*ymed(j)/(G*MSTAR);
	lapler = ((e[llyp]/rho[llyp]-e[ll]/rho[ll])*SurfY(j+1,k)-	\
		(e[ll]/rho[ll]-e[llym]/rho[llym])*SurfY(j,k))/(ymin(j+1)-ymin(j));
	lapler += ((e[llxp]/rho[llxp]-e[ll]/rho[ll])*SurfX(j,k)-	\
		   (e[ll]/rho[ll]-e[llxm]/rho[llxm])*SurfX(j,k))/zone_size_x(j,k);
	lapler *= InvVol(j,k)*chi;
	lapl[ll] = lapler;
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
