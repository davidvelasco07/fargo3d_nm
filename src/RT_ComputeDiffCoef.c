//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ComputeDiffusionCoefficients_cpu () {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Energyrad);
  INPUT(OpaPlanck);
  OUTPUT(DiffCoef);
//<\USER_DEFINED>


//<EXTERNAL>
  real* dens     = Density->field_cpu;
  real* enrad    = Energyrad->field_cpu;
  real* opar     = OpaPlanck->field_cpu;
  real* diffcoef = DiffCoef->field_cpu;
  int pitch      = Pitch_cpu;
  int stride     = Stride_cpu;
  int size_x     = XIP;
  int size_y     = Ny+2*NGHY-1;
  int size_z     = Nz+2*NGHZ-1;
  real dx        = Dx;
//<\EXTERNAL>

//<INTERNAL>
  real graderrad;
  real graderazi;
  real gradercol;
  real erre;
  real lambda;
  real coef;
  int i;
  int j;
  int k;
  int ll;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
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
      for (i=XIM; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
	coef = 1.0/(opar[ll]*dens[ll]);
	graderazi = .5*(enrad[lxp]-enrad[lxm])/zone_size_x(j,k);
	graderrad = .5*(enrad[lyp]-enrad[lym])/zone_size_y(j,k);
	gradercol = .5*(enrad[lzp]-enrad[lzm])/zone_size_z(j,k);
	erre = sqrt(graderrad*graderrad+graderazi*graderazi+gradercol*gradercol) \
	  *coef/enrad[ll];

	// Flux limiter implemented below
	/////////////
	if (erre <=2.) {
	  lambda = 2./(3.+sqrt(9. +10.*erre*erre));
	} else {
	  lambda = 10./(10.*erre+9.+sqrt(81.+180.*erre));
	}
	// End of flux limiter
	/////////////
	diffcoef[ll] = lambda*coef*C0;
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
//  comm_cpu (DIFFCOEF);  This should be needed when we add the source terms
//<\LAST_BLOCK>
}
