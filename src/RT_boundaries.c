//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void SetRTBoundaryConditions_cpu () {

//<USER_DEFINED>
  INPUT(LeftY);
  INPUT(RightY);
  INPUT(RightZ);
  INPUT(BB);
  INPUT(Energyrad);
  INPUT(DiffCoef);
  OUTPUT(LeftY);
  OUTPUT(RightY);
  OUTPUT(RightZ);
  OUTPUT(BB);
  OUTPUT(Energyrad);
  OUTPUT(DiffCoef);
//<\USER_DEFINED>


//<EXTERNAL>
  real* lefty  = LeftY->field_cpu;
  real* righty = RightY->field_cpu;
  real* rightz = RightZ->field_cpu;
  real* bb     = BB->field_cpu;
  real* erad   = Energyrad->field_cpu;
  real* diffcoef = DiffCoef->field_cpu;
  int pitch    = Pitch_cpu;
  int stride   = Stride_cpu;
  int size_x   = Nx+2*NGHX;
  int size_y   = Ny+2*NGHY;
  int size_z   = Nz+2*NGHZ;
  real cpuymin = Ymin[NGHY];
  real cpuzmin = Zmin[NGHZ];
  real ymin    = corner_min0[1];
  real zmin    = corner_min0[2];
  real cpuymax = Ymin[Ny+NGHY];
  real cpuzmax = Zmin[Nz+NGHZ];
  real ymax    = corner_max0[1];
  real zmax    = corner_max0[2];
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
//<\INTERNAL>
  
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

	// Boundary condition at the inner radius (equivalent to
	// putting a mirror there that reflects outgoing photons)
	if ((cpuymin == ymin) && (j == NGHY)) {
	  bb[l] += lefty[l];
	  lefty[l] = 0.0;
	}

	// Boundary condition at the outer radius (equivalent to
	// putting a mirror there that reflects outgoing photons)
	if ((cpuymax == ymax) && (j == size_y-NGHY-1)) {
	  bb[l] += righty[l];
	  righty[l] = 0.0;
	}

	// Boundary condition at the equator for a half-disk (again
	// equivalent to having a mirror).
#ifdef HALFDISK
	if ((cpuzmax == zmax) && (k == size_z-NGHZ-1)) {
	  bb[l] += rightz[l];
	  rightz[l] = 0.0;
	}
#endif

	if ((cpuzmin == zmin) && (k == 0)) {
	  diffcoef[l] = diffcoef[lzp];
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
