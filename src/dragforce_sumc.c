//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void DragForce_SumC_cpu (real dt) {
 
//<USER_DEFINED>
  INPUT(Density);
  INPUT(Qs);
  INPUT(DensStar);
  OUTPUT(DensStar);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* c    = DensStar->field_cpu;
  real* pref = Qs->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real invstokesnumber = Coeffval[0];
  real invparticlesize = Coeffval[1];
  real rhosolid        = Coeffval[2];
  real tslim           = TSLIM;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real alphak;
  real sk;
//<\INTERNAL>

//<CONSTANT>
// real ymin(Ny+2*NGHY+1);
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
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;

#ifdef STOKESNUMBER
	alphak  = pref[ll]*invstokesnumber;
#endif
#ifdef DUSTSIZE
	alphak  = max2( pref[ll]*sqrt(8./M_PI)*invparticlesize/rhosolid, sqrt(G*MSTAR/(ymed(j)*ymed(j)*ymed(j)))/tslim )  ;
#endif
	sk      = dt*alphak/(1+dt*alphak);

	if (fluidtype == GAS)  sk = 1.0;
	c[ll] += dens[ll]*sk;

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
