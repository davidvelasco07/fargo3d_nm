//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
#include "J_jupiter.h"
//<\INCLUDES>

void MatrixIterate_cpu (int parity, real omega) {

//<USER_DEFINED>
  INPUT (LeftX);
  INPUT (LeftY);
  INPUT (LeftZ);
  INPUT (RightX);
  INPUT (RightY);
  INPUT (RightZ);
  INPUT (BB);
  INPUT (RHS);
  INPUT (Energyrad);
  OUTPUT (Resid);
  OUTPUT (Energyrad);
//<\USER_DEFINED>


//<EXTERNAL>
  real* leftx  = LeftX->field_cpu;
  real* lefty  = LeftY->field_cpu;
  real* leftz  = LeftZ->field_cpu;
  real* rightx = RightX->field_cpu;
  real* righty = RightY->field_cpu;
  real* rightz = RightZ->field_cpu;
  real* bb     = BB->field_cpu;
  real* rhs    = RHS->field_cpu;
  real* resid  = Resid->field_cpu;
  real* enrad  = Energyrad->field_cpu;
  char* hidden = Current_Jupiter_Patch->Hidden;
  int pitch    = Pitch_cpu;
  int stride   = Stride_cpu;
  int size_x   = Nx+NGHX;
  int size_y   = Ny+NGHY;
  int size_z   = Nz+NGHZ;
  int x0       = Current_Jupiter_Patch->pcorner_min[0];
  int y0       = Current_Jupiter_Patch->pcorner_min[1];
  int z0       = Current_Jupiter_Patch->pcorner_min[2];
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real res;
//<\INTERNAL>
  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=NGHZ; k<size_z; k++) {
#endif
#ifdef Y
    for (j=NGHY; j<size_y; j++) {
#endif
#ifdef X
      for (i=NGHX; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
	if ((i+j+k+x0+y0+z0) % 2 == parity && !hidden[ll]) { //Checkerboard avoids race condition
	  res = bb[ll]*enrad[ll]-rhs[ll] +
	    leftx[ll]*enrad[lxm] + rightx[ll]*enrad[lxp] +
	    lefty[ll]*enrad[lym] + righty[ll]*enrad[lyp] +
	    leftz[ll]*enrad[lzm] + rightz[ll]*enrad[lzp];	
	  enrad[ll] -= omega*res/bb[ll];
	  resid[ll] = fabs(res);
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
