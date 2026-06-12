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
  real constant = (1.-GAMMA)*THDIFFUSIVITY/R_MU;
  real dx = Dx;
  real rho0 = RHO0;
  real e0 = CS*CS*RHO0/(GAMMA*(GAMMA-1.0));
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
  real e_a;
  real rho_a;
  real chi_a;
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
	

#if THDIFFUSION == 0 //Constant thermal diffusivity
#ifdef X
	llxm = lxm;
	fx[ll] = constant * SurfX(j,k) * (e_rho - e[llxm]/rho[llxm]) / zone_size_x(j,k);
#endif
#ifdef Y
	llym = lym;
	fy[ll] = constant * SurfY(j,k) * (e_rho - e[llym]/rho[llym]) / zone_size_y(j,k);
#endif
#ifdef Z
	llzm = lzm;
	fz[ll] = constant * SurfZ(j,k) * (e_rho - e[llzm]/rho[llzm]) / zone_size_z(j,k);
#endif
#endif
#if THDIFFUSION == 1 //Chi scales as the inverse of the density
#ifdef X
	llxm = lxm;
	rho_a = 2*rho0/(rho[ll]+rho[llxm]);
	fx[ll] = constant * SurfX(j,k) * rho_a * (e_rho - e[llxm]/rho[llxm]) / zone_size_x(j,k);
#endif
#ifdef Y
	llym = lym;
	rho_a = 2*rho0/(rho[ll]+rho[llym]);
	fy[ll] = constant * SurfY(j,k) * rho_a * (e_rho - e[llym]/rho[llym]) / zone_size_y(j,k);
#endif
#ifdef Z
	llzm = lzm;
	rho_a = 2*rho0/(rho[ll]+rho[llzm]);
	fz[ll] = constant * SurfZ(j,k) * rho_a * (e_rho - e[llzm]/rho[llzm]) / zone_size_z(j,k);
#endif
#endif
#if THDIFFUSION == 2 //Chi scales as T^3/rho^2
#ifdef X
	llxm = lxm;
       	e_a = 0.5*(e[ll]+e[llxm])/e0;
	rho_a = 2*rho0/(rho[ll]+rho[llxm]);
	chi_a = e_a*e_a*e_a*rho_a*rho_a*rho_a*rho_a*rho_a;
	fx[ll] = constant * SurfX(j,k) * chi_a * (e_rho - e[llxm]/rho[llxm]) / zone_size_x(j,k);
#endif
#ifdef Y
	llym = lym;
	e_a = 0.5*(e[ll]+e[llym])/e0;
	rho_a = 2*rho0/(rho[ll]+rho[llym]);
	chi_a = e_a*e_a*e_a*rho_a*rho_a*rho_a*rho_a*rho_a;
	fy[ll] = constant * SurfY(j,k) * chi_a * (e_rho - e[llym]/rho[llym]) / zone_size_y(j,k);
#endif
#ifdef Z
	llzm = lzm;
	e_a = 0.5*(e[ll]+e[llzm])/e0;
	rho_a = 2*rho0/(rho[ll]+rho[llzm]);
	chi_a = e_a*e_a*e_a*rho_a*rho_a*rho_a*rho_a*rho_a;
	fz[ll] = constant * SurfZ(j,k) * chi_a * (e_rho - e[llzm]/rho[llzm]) / zone_size_z(j,k);
#endif
#endif
#if THDIFFUSION == 3 //Planet-disc:Constant thermal diffusivity
#ifdef X
	llxm = lxm;
	fx[ll] = constant * SurfX(j,k) * 0.5 * (rho[ll]+rho[llxm]) * (e_rho - e[llxm]/rho[llxm]) / zone_size_x(j,k);
#endif
#ifdef Y
	llym = lym;
	fy[ll] = constant * SurfY(j,k) * 0.5 * (rho[ll]+rho[llym]) * (e_rho - e[llym]/rho[llym]) / zone_size_y(j,k);
#endif
#ifdef Z
	llzm = lzm;
	fz[ll] = constant * SurfZ(j,k) * 0.5 * (rho[ll]+rho[llzm]) *  (e_rho - e[llzm]/rho[llzm]) / zone_size_z(j,k);
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
