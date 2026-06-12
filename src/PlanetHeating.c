//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void PlanetHeating(real dt) {
  int n;
  real heat= dt*HEATINGRATE;
  real* xp = Sys->x_cpu;
  real* yp = Sys->y_cpu;
  real* zp = Sys->z_cpu;
  int nb = Sys->nb;
  real xs,ys,zs,dx,dy,dz;
  dx=(Xmed(1)-Xmed(0));
  dy=(Ymed(1)-Ymed(0));
  dz=(Zmed(1)-Zmed(0));
  
  for (n=0;n<nb;n++){
#ifdef PLANET_HEATING
#if PLANET_HEATING==0
    xs = atan2(yp[n],xp[n]);
    ys = sqrt(xp[n]*xp[n]+yp[n]*yp[n]+zp[n]*zp[n]);
    zs = atan2(sqrt(xp[n]*xp[n]+yp[n]*yp[n]),zp[n]);
    
    FARGO_SAFE(_PlanetHeating(xs,ys,zs,dx,dy,dz,heat));
#endif
#if PLANET_HEATING==1
    FARGO_SAFE(_PlanetHeating(xp[n],yp[n],zp[n],dx,dy,dz,heat));
#endif
#endif
  }
}

void _PlanetHeating_cpu(real xp, real yp, real zp, real dx, real dy, real dz, real heat) {

//<USER_DEFINED>
  INPUT(Energy);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* e   = Energy->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX-1;
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
  int i0 = 0;
  int j0 = 0;
  int k0 = 0;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int n;
  real x;
  real y;
  real z;
  real frac;
  real r;
  real planetdistance;
  real smoothing;
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
// real ASPECTRATIO(1);
// real FLARINGINDEX(1);
// real THICKNESSSMOOTHING(1);
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
#ifdef PLANET_HEATING
#if PLANET_HEATING==0
	  frac=1;
	  
	  x=fabs(xp-xmed(i));
	  y=fabs(yp-ymed(j));
	  z=fabs(zp-zmed(k));
	  if(x<dx&& y<dy && z<dz){
	    frac *= 1-x/dx;
	    frac *= 1-y/dy;
	    frac *= 1-z/dz;
	    e[l] += heat*InvVol(j, k)*frac;
	  }
#endif
#if PLANET_HEATING==1
	  planetdistance = sqrt(xp*xp+yp*yp+zp*zp);
	  smoothing = ASPECTRATIO* pow(planetdistance/R0,FLARINGINDEX)*
	    planetdistance*THICKNESSSMOOTHING;
	  
	  x=fabs(xp-XC);
	  y=fabs(yp-YC);
	  z=fabs(zp-ZC);
	  r = sqrt(x*x+y*y*+z*z);
	  
	  e[l] += heat*exp(-r/smoothing)/(8*M_PI*smoothing*smoothing*smoothing);
#endif
#endif
	  //}
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

