//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ReMap_cpu(real inner, real outer) { 

  /*inner = new xmin, outer = new xmax   */

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Vx);
  INPUT(Vy);
  INPUT(Energy);
  OUTPUT(Pressure);
  OUTPUT(Mpx);
  OUTPUT(Mmx);
  OUTPUT(Mmy);
//<\USER_DEFINED>

//<EXTERNAL>
  int size_x = Nx+2*NGHX;
  int size_y= Ny+2*NGHY;
  int size_z= Nz+2*NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real* rho = Density->field_cpu;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
  real* cs = Energy->field_cpu;
  real* rho_temp = Pressure->field_cpu;
  real* vx_temp = Mpx->field_cpu;
  real* vy_temp = Mmx->field_cpu;
  real* cs_temp = Mmy->field_cpu;
  real ny = NY;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real alpha;
  real r;
  real omega;
  int jdown;
  int jdown_med;
  real ynew;
  real ynew_med;
  real fact; 
  real deno;
  real delta_y = (outer-inner)/ny;
//<\INTERNAL>

//<CONSTANT>
// real ymin(Ny+2*NGHY+1);
// real FLARINGINDEX(1);
// real OMEGAFRAME(1);
// real SIGMASLOPE(1);
// real SIGMA0(1);
// real Y0(1);
// real ASPECTRATIO(1);
//<\CONSTANT>

  i = j = k = 0;

//<MAIN_LOOP>
#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif

//<#>

	ll = l;

	ynew_med = 0.5*(inner + delta_y * (j+Y0-NGHY) + inner + delta_y * (j+1+Y0-NGHY));
	jdown_med = (int)((ynew_med-ymed(0))/(ymed(size_y-1)-ymed(0))*(size_y-1));

	ynew = inner + delta_y * (j+Y0-NGHY);	
	jdown = (int)((ynew-ymin(0))/(ymin(size_y)-ymin(0))*(size_y));

	//Staggered fields	
	if (jdown < size_y-1 && jdown >= 0) { //Indices must keep ghosts! (MPI)

	  fact = (ynew-ymin(jdown));
	  deno = 1.0/(ymin(jdown+1)-ymin(jdown));

	  alpha = (vy[i+(jdown+1)*pitch+k*stride]-vy[i+(jdown)*pitch+k*stride])*deno;
	  vy_temp[ll] = alpha*fact + vy[i+(jdown)*pitch+k*stride];
	}
	else {
	  vy_temp[ll]  = 0.0;
	}
	
	//Centered fields
	if (jdown_med < size_y-1 && jdown_med >= 0) { //Indices must keep ghosts! (MPI)

	  fact = (ynew_med-ymed(jdown_med));
	  deno = 1.0/(ymed(jdown_med+1)-ymed(jdown_med));

	  alpha = (rho[i+(jdown_med+1)*pitch+k*stride]-rho[i+(jdown_med)*pitch+k*stride])*deno;
	  rho_temp[ll] = alpha*fact + rho[i+(jdown_med)*pitch+k*stride];
	  
	  alpha = (vx[i+(jdown_med+1)*pitch+k*stride]-vx[i+(jdown_med)*pitch+k*stride])*deno;
	  vx_temp[ll] = alpha*fact + vx[i+(jdown_med)*pitch+k*stride];

	  alpha = (cs[i+(jdown_med+1)*pitch+k*stride]-cs[i+(jdown_med)*pitch+k*stride])*deno;
	  cs_temp[ll] = alpha*fact + cs[i+(jdown_med)*pitch+k*stride];
	}
	else {
	  r = ynew_med;
	  omega  = sqrt(G*MSTAR/r/r/r);
	  rho_temp[ll] = SIGMA0*pow((r/R0),-SIGMASLOPE);
	  cs_temp[ll]  = sqrt(G*MSTAR/r)*ASPECTRATIO*pow(r/R0, FLARINGINDEX);
	  vx_temp[ll]  = omega*r*sqrt(1.0+pow(ASPECTRATIO,2) * 
				     pow(r/R0,2*FLARINGINDEX) * 
				     (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE)) - OMEGAFRAME*r;
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

//<LAST_BLOCK>
  
  ReMapGeometryStock(inner, outer);
  ReMapCopy();

//<\LAST_BLOCK>

}
