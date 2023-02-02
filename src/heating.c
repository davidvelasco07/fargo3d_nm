//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void compute_planetheating(real dt) {
  int n;
  int i;
  real heat;
  real* xp = Sys->x;
  real* yp = Sys->y;
  real* zp = Sys->z;
  real* Mp = Sys->mass;
  real rhop = 706118396.1461967; //3g/cm^3 -> c.u
  real Rp;
  int nb = Sys->nb;  
  for (n=0;n<nb;n++){
    heat=0;
    for(i=0;i<NbFluids;i++){
        //L = G Mp M'p/Rp
        Rp = pow(.25*3*Mp[n]/(M_PI*rhop),1./3);
        heat += G*Mp[n]*M_dot[FluidIndex][n]/Rp;
    }
    FARGO_SAFE(PlanetHeating(xp[n],yp[n],zp[n],heat*dt));
  }
}

void PlanetHeating_cpu(real xp, real yp, real zp, real heat) {

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
  real dx;
  real dy;
  real dz;
  real dxp;
  real dyp;
  real dzp;
  real frac;
  real dist;
  real rp;
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
#if PLANET_HEATING==0
	      if(dxp<dx && dyp<dy && dzp<dz){
	        frac *= 1-dxp/dx;
	        frac *= 1-dyp/dy;
	        frac *= 1-dzp/dz;
	        e[l] += heat*InvVol(j, k)*frac;
	      } 
#endif
#if PLANET_HEATING==1
        rp = sqrt(xp*xp+yp*yp);
        smoothing = ASPECTRATIO*pow(rp/R0,FLARINGINDEX)*rp*THICKNESSSMOOTHING;
        dist = sqrt(dxp*dxp+dyp*dyp+dzp*dzp);
        frac = exp(-dist/smoothing)/(8*M_PI*smoothing*smoothing*smoothing);
	      e[l] += heat*frac*InvVol(j, k);
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

