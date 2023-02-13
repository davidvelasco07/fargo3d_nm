//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void SubStep4_a_cpu (real dt) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Energy);
  INPUT(ThermalDiff);
#ifdef X
  OUTPUT(Mmx);
#endif
#ifdef Y
  OUTPUT(Mmy);
#endif
#ifdef Z
  OUTPUT(Mmz);
#endif

//<\USER_DEFINED>

//<EXTERNAL>
  real* rho = Density->field_cpu;
  real* e = Energy->field_cpu;
  real* thd = ThermalDiff->field_cpu;
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
  int size_x = Nx+2*NGHX; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real constant = (1.-GAMMA)/R_MU;
  real dx = Dx;
//<\EXTERNAL>
  
//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
#ifdef X
  int llxm;
#endif
#ifdef Y
  int llym;
#endif
#ifdef Z
  int llzm;
#endif
  real e_rho;
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
	e_rho = e[ll]/rho[ll];
#ifdef THDIFFUSION
#ifdef X
	llxm = lxm;
	fx[ll] = constant * 0.25*(thd[ll]+thd[llxm])*(rho[ll]+rho[llxm])*(e_rho-e[llxm]/rho[llxm]) * SurfX(j,k)/zone_size_x(j,k);
#endif
#ifdef Y
	llym = lym;
	fy[ll] = constant * 0.25*(thd[ll]+thd[llym])*(rho[ll]+rho[llym])*(e_rho-e[llym]/rho[llym]) * SurfY(j,k)/zone_size_y(j,k);
#endif
#ifdef Z
	llzm = lzm;
	fz[ll] = constant * 0.25*(thd[ll]+thd[llzm])*(rho[ll]+rho[llzm])*(e_rho-e[llzm]/rho[llzm]) * SurfZ(j,k)/zone_size_z(j,k);
#endif
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
