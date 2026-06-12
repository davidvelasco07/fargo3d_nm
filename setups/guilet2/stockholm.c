//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void StockholmBoundary_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT2D(Density0);
  OUTPUT(Density);
#ifdef ADIABATIC
  INPUT(Energy);
  INPUT2D(Energy0);
  OUTPUT(Energy);
#endif
#ifdef X
  INPUT(Vx);
  INPUT2D(Vx0);
  OUTPUT(Vx);
#endif
#ifdef Y
  INPUT(Vy);
  INPUT2D(Vy0);
  OUTPUT(Vy);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho  = Density->field_cpu;
  real* rho0 = Density0->field_cpu;
#ifdef X
  real* vx  = Vx->field_cpu;
  real* vx0 = Vx0->field_cpu;
#endif
#ifdef Y
  real* vy  = Vy->field_cpu;
  real* vy0 = Vy0->field_cpu;
#endif
#ifdef ADIABATIC
  real* e    = Energy->field_cpu;
  real* e0   = Energy0->field_cpu;
#endif
  int pitch   = Pitch_cpu;
  int stride  = Stride_cpu;
  int size_x  = Nx;
  int size_y  = Ny+2*NGHY;
  int size_z  = Nz+2*NGHZ;
  int pitch2d = Pitch2D;
  real y_min = YMIN;
  real y_max = YMAX;
  real r0 = R0;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  // Below: we match the standard Stockholm's prescription
  // see De Val Borro et al. (2006), section 3.2
  //  real Y_inf = y_min + (y_max-y_min)*0.0476;
  //  real Y_sup = y_max - (y_max-y_min)*0.19;
  real Y_inf = y_min + (y_max-y_min)*0.05;
  real Y_sup = y_max - (y_max-y_min)*0.10;
  real rampy;
  real tau;
  real ds = 0.1;
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
	rampy = 0.0;
#ifdef Y
	if(ymed(j) > Y_sup) {
	  rampy   = (ymed(j)-Y_sup)/(y_max-Y_sup);
	}
	if(ymed(j) < Y_inf) {
	  rampy   = (Y_inf-ymed(j))/(Y_inf-y_min);
	}
	rampy *= rampy;		/* Parabolic ramp as in De Val Borro et al (2006) */
#endif
	tau = ds*sqrt(ymed(j)*ymed(j)*ymed(j)/G/MSTAR);
	if(rampy>0.0) {
	  tau = tau/rampy;
	  //	  rho[l] = (rho[l]*tau+rho0[l2D]*dt)/(dt+tau);
#ifdef X
	  vx[l] = (vx[l]*tau+vx0[l2D]*dt)/(dt+tau);
#endif
#ifdef Y
	  vy[l] = (vy[l]*tau+vy0[l2D]*dt)/(dt+tau);
#endif
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
