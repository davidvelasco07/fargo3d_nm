//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void DragForce_Coeff_cpu () {

  /*This function reutns the prefactor of the drag coeffcient which
    is independent of the fluid type. The Stokes number or particle size
    and solid density are added through alpha[n] in each source file */
  
//<USER_DEFINED>
#ifdef DUSTSIZE
  INPUT(QR);
  INPUT(Energy);
#endif
  OUTPUT(Qs);
//<\USER_DEFINED>

//<EXTERNAL>
  real* coeff      = Qs->field_cpu;
#ifdef DUSTSIZE
  real* dens_gas   = QR->field_cpu;
  real* energy_gas = Energy->field_cpu; 
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real omega;
  real cs;
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
	omega    = sqrt(G*MSTAR/(ymed(j)*ymed(j)*ymed(j)));
	coeff[ll] = omega;
#endif
	
#ifdef DUSTSIZE
#ifdef ADIABATIC
	cs = sqrt(GAMMA*(GAMMA-1)*energy_gas[ll]/dens_gas[ll]);
#endif
#ifdef ISOTHERMAL
	cs = energy_gas[ll];
#endif
	coeff[ll] = dens_gas[ll]*cs;
#endif
	
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
