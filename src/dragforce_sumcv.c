//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void DragForce_SumCV(real dt, int option) {

  if (option == 0)
    _DragForce_SumCV(dt, 1, 0, 0, Vx_temp, Mmx);
  if (option == 1)
    _DragForce_SumCV(dt, 0, 1, 0, Vy_temp, Mmy);
  if (option == 2)
    _DragForce_SumCV(dt, 0, 0, 1, Vz_temp, Mmz);
  
}

void _DragForce_SumCV_cpu(real dt, int idx, int idy, int idz, Field *V, Field *Cv) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Qs);
  INPUT(V);
  INPUT(Cv);
  OUTPUT(Cv);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* pref = Qs->field_cpu;
  real* v    = V->field_cpu;
  real* cv   = Cv->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real invstokesnumber = Coeffval[0];
  real invparticlesize = Coeffval[1];
  real rhosolid        = Coeffval[2];
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  int lm;
  real alphak;
  real sk;
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
	
#ifdef STOKESNUMBER
	alphak  = 0.5*(pref[ll]+pref[lm])*invstokesnumber;
#endif
#ifdef DUSTSIZE
	alphak  = 0.5*(pref[ll]+pref[lm])*sqrt(8./M_PI)*invparticlesize/rhosolid;
#endif
	sk      = dt*alphak/(1+dt*alphak);

	if (fluidtype == GAS)  sk = 1.0;

	cv[ll] += 0.5*(dens[ll]+dens[lm])*v[ll]*sk;

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
