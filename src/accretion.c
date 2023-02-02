//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void compute_accretion(real dt) {
  int i,j,k,n;
  real* xp = Sys->x;
  real* yp = Sys->y;
  real* zp = Sys->z;
  real* vxp = Sys->vx;
  real* vyp = Sys->vy;
  real* vzp = Sys->vz;
  real* mp = Sys->mass;
  int nb = Sys->nb;

  for (n=0;n<nb;n++){
    Accretion(xp[n], yp[n], zp[n]);
    M_dot[FluidIndex][n] = reduction_full_SUM(Energy, NGHY, Ny+NGHY, 0, Nz+2*NGHZ);
    M_acc[FluidIndex][n] += M_dot[FluidIndex][n];
  }
}

void Accretion_cpu (real xp, real yp, real zp) {
//<USER_DEFINED>
  INPUT(Density);
  OUTPUT(Density);
  OUTPUT(Energy);
#ifdef X
  INPUT(Vx);
#endif
#ifdef Y
  INPUT(Vy);
#endif
#ifdef Z
  INPUT(Vz);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho_d = Density->field_cpu;
#ifdef X 
  real* vx_d = Vx->field_cpu;
#endif
#ifdef Y
  real* vy_d = Vy->field_cpu;
#endif
#ifdef Z
  real* vz_d = Vz->field_cpu;
#endif
  real* accreted = Energy->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real stokes = 1./Coeffval[0];
//<\EXTERNAL>
  
//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
  int ip;
  int jp;
  int kp;
  int lp;
  real x;
  real y;
  real z;
  real vx;
  real vy;
  real vz;
  real dx;
  real dy;
  real dz;
  real dxp;
  real dyp;
  real dzp;
  real r_polar;
  real dist;
  real frac;
  real removed;
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
// real OMEGAFRAME(1);
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
  r_polar = sqrt(xp*xp+yp*yp);
#ifdef CARTESIAN
  x = xp;
  y = yp;
  z = zp;
#endif
#ifdef CYLINDRICAL
  x = atan2(yp,xp);
  y = sqrt(xp*xp+yp*yp);
  z = zp;
#endif
#ifdef SPHERICAL
  x = atan2(yp,xp);
  y = sqrt(xp*xp+yp*yp+zp*zp);
  z = atan2(sqrt(xp*xp+yp*yp),zp);
#endif
  dx = xmed(1)-xmed(0);
  dy = ymed(1)-ymed(0);
  dz = zmed(1)-zmed(0);
  dxp=fabs(x-xmed(i));
	dyp=fabs(y-ymed(j));
	dzp=fabs(z-zmed(k));
  frac=1;
  //Just in case let's clean whatever was in this field
  accreted[ll] = 0;
	if(dxp<dx && dyp<dy && dzp<dz){
	    frac *= 1-dxp/dx;
	    frac *= 1-dyp/dy;
	    frac *= 1-dzp/dz;
      removed = frac*rho_d[ll];
      rho_d[ll] -= removed;
      accreted[ll] = Vol(j,k)*removed;
      //printf("(%d,%d,%d)->%g, dens=%lf, removed=%lf, Vol=%.18g, Macc=%.18g\n",i,j,k,frac,rho_d[ll],removed,Vol(j,k),accreted[ll]);
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
