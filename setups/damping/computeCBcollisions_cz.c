//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ComputeCBcollisions_cz_cpu (real dt) {
 
//<USER_DEFINED>
  INPUT(Density);
  INPUT(Total_Density);
  INPUT(Qs);
  INPUT(Vz_temp);
  OUTPUT(Mpz);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens     = Density->field_cpu;
  real* dens_gas = Total_Density->field_cpu;
  real* pref     = Qs->field_cpu;
  real* vz       = Vz_temp->field_cpu;
  real* cz       = Mpz->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int n     = FluidIndex;
  int fluidtype = Fluidtype;
  real* alpha = Alpha;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real gamma_k;
  real s_k;
  real epsilon;
  real _cz;
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

	epsilon = (dens[ll]+dens[lzm])/(dens_gas[ll]+dens_gas[lzm]);

	gamma_k = 0.5*(pref[ll]+pref[lzm]);
	s_k     = dt*gamma_k/(1+dt*gamma_k);

	if (fluidtype == GAS)  _cz = vz[ll];
	else _cz = s_k*epsilon*vz[ll];
	cz[ll] += _cz;
	
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
