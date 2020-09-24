//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ComputeCBcollisions_cv(real dt, int option) {

  if (option == 0)
    _ComputeCBcollisions_cv(dt, 1, 0, 0, Vx_temp, Mpx);
  if (option == 1)
    _ComputeCBcollisions_cv(dt, 0, 1, 0, Vy_temp, Mpy);
  if (option == 2)
    _ComputeCBcollisions_cv(dt, 0, 0, 1, Vz_temp, Mpz);
  
}

void _ComputeCBcollisions_cv_cpu(real dt, int idx, int idy, int idz, Field *V, Field *Cv) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Total_Density);
  INPUT(Qs);
  INPUT(V);
  OUTPUT(Cv);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens     = Density->field_cpu;
  real* dens_gas = Total_Density->field_cpu;
  real* pref     = Qs->field_cpu;
  real* v        = V->field_cpu;
  real* cv       = Cv->field_cpu;
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
  real _cv;
  real epsilon;
//<\INTERNAL>

//<CONSTANT>
// real Alpha(NFLUIDS*NFLUIDS);
//<\CONSTANT>

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
	
	epsilon = (dens[ll]+dens[lm])/(dens_gas[ll]+dens_gas[lm]);
	gamma_k = 0.5*(pref[ll]+pref[lm]);
	s_k     = dt*gamma_k/(1+dt*gamma_k);

	if (fluidtype == GAS)  _cv = v[ll];
	else _cv  = s_k*epsilon*v[ll];
	
	cv[ll] += _cv;
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
