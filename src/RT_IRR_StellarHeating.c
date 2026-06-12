//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ComputeStellarHeating_cpu () {

//<USER_DEFINED>
  INPUT(Tau);
  OUTPUT(StellarRad);

  /* Under the standard practice of FARGO3D, RSTAR and TSTAR should
     have scaling rules and should be expressed in the current unit
     system in a given parameter file. However, this would be hardly
     legible as we use dimensionless variables for most of the
     parameters, which are therefore O(1). We prefer here to enforce
     their definition in cgs, and we convert in this function on the
     fly. Harmless warnings are issued at build time. */
  real rstar = RSTAR/R0_CGS*R0;
  real tstar = TSTAR*(G*MSTAR/R0/R_MU)/(G_CGS*MSTAR_CGS/R0_CGS/R_MU_CGS);
//<\USER_DEFINED>

//<EXTERNAL>
  real* tau        = Tau->field_cpu;
  real* stellarrad = StellarRad->field_cpu;
  real  fstar      = STEFANK * rstar*rstar * tstar*tstar*tstar*tstar;
  int pitch        = Pitch_cpu;
  int stride       = Stride_cpu;
  int size_x       = Nx+2*NGHX;
  int size_y       = Ny+2*NGHY;
  int size_z       = Nz+2*NGHZ;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
//<\INTERNAL>

//<CONSTANT>
// real Syj(Ny+2*NGHY);
// real Syk(Nz+2*NGHZ);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real InvVj(Ny+2*NGHY);
// real HEATTAPERING(1);
// real HEATSHIELD(1);
// real HEATRADIUS(1);
//<\CONSTANT>

  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=1; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
	if ((ymed(j) > HEATRADIUS) && ((zmed(k) <= .5*M_PI-HEATSHIELD) || (zmed(k) >= .5*M_PI+HEATSHIELD))) { 	  //Non-divergent test
	  stellarrad[ll] = fstar*SurfY(j,k)/(Vol(j,k)*ymin(j)*ymin(j))*(exp(-tau[lym])-exp(-tau[ll]));
	  if (ymed(j) < HEATRADIUS+HEATTAPERING)
	    stellarrad[ll] *= pow(sin(.5*M_PI*((ymed(j)-HEATRADIUS)/HEATTAPERING)),2.0);
	}
	else
	  stellarrad[ll] = 0.0;
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
