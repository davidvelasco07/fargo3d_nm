//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void UpdateVelcollisions(real dt, int option) {

  if (option == 0)
    _UpdateVelcollisions(dt, 1, 0, 0, Vx_temp, Mpx);
  if (option == 1)
    _UpdateVelcollisions(dt, 0, 1, 0, Vy_temp, Mpy);
  if (option == 2)
    _UpdateVelcollisions(dt, 0, 0, 1, Vz_temp, Mpz);
  
}

void _UpdateVelcollisions_cpu(real dt, int idx, int idy, int idz, Field *V, Field *Cv) {

//<USER_DEFINED>
  INPUT(Slope);
  INPUT(Cv);
  INPUT(Qs);
  INPUT(V);
  OUTPUT(V);
//<\USER_DEFINED>

//<EXTERNAL>
  real* v    = V->field_cpu;
  real* cv   = Cv->field_cpu;
  real* c    = Slope->field_cpu;
  real* pref = Qs->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  int lm;
  real gamma_k;
  real s_k;
  real dst;
//<\INTERNAL>
  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=1; k<size_z; k++) {
#endif
#ifdef Y
    for (j=1; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;

	lm = idx*lxm + idy*lym + idz*lzm;
	
	gamma_k = 0.5*(pref[ll]+pref[lm]);
	s_k     = dt*gamma_k/(1+dt*gamma_k);
	
	if (fluidtype == GAS)  v[ll] = cv[ll]/( 1. + 0.5*(c[ll]+c[lm]) );
	else                   v[ll] = s_k*cv[ll]/(1. + 0.5*(c[ll]+c[lm]) ) + v[ll]/(1.+ dt*gamma_k);
	
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
