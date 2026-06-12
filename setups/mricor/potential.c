//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void compute_potential(real dt) {

  if (strcmp (PLANETCONFIG, "NONE") != 0) {
    mastererr ("This setup is not meant to include an evolving planetary system");
    mastererr ("A planet with mass PLANETMASS is set on a fixed circular orbit at R0");
    exit (1);
  }
  FARGO_SAFE(Potential ());
}

void Potential_cpu() {
  
//<USER_DEFINED>
  OUTPUT(Pot);
  real planetmass_taper;
  if (MASSTAPER == 0.0)
    planetmass_taper = 1.0;
  else
    planetmass_taper = (PhysicalTime - PhysicalTimeRestart >= MASSTAPER ? 1.0 : .5*(1.0-cos(M_PI*(PhysicalTime-PhysicalTimeRestart)/MASSTAPER)));
//<\USER_DEFINED>

//<EXTERNAL>
  real* pot  = Pot->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real smoothing = 3.0*R0*THICKNESSSMOOTHING*ASPECTRATIO;
  real planetmass = PLANETMASS*planetmass_taper;
//<\EXTERNAL>
  
//<INTERNAL>
  int i;
  int j;
  int k;
  real dist;
  real phi;
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
	phi =  -G*MSTAR/ymed(j); //Potential from star

	dist = (YC-3.0*R0)*(YC-3.0*R0)+XC*XC; // Planet at 90 degrees owing to limits
	pot[l] = phi-G*planetmass/sqrt(dist+smoothing*smoothing); //Potential from planets
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
