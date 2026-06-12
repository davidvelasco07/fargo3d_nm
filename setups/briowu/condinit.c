#include "fargo3d.h"

void CondInit(){
  
  int i,j,k;
  
  real* bx = Bx->field_cpu;
  real* by = By->field_cpu;
  real* bz = Bz->field_cpu;
  real* rho = Density->field_cpu;
  real* cs  = Energy->field_cpu;
  real* vy = Vy->field_cpu;
  real* vx = Vx->field_cpu;
  real* vz = Vz->field_cpu;

  i = j = k = 0;

#ifdef Z
  for (k=0; k<Nz+2*NGHZ; k++) {
#endif
#ifdef Y
    for (j=0; j<Ny+2*NGHY; j++) {
#endif
#ifdef X
      for (i=0; i<Nx; i++) {
#endif
	cs[l]  = 1.0/(GAMMA-1.0);
	rho[l] = 1.0;
	bz[l]  = 0.75;
	by[l]  = 1.0;
	if(Zmed(k)>0.5*(ZMAX+ZMIN)) {
	  cs[l]  = 0.1/(GAMMA-1.0);
	  rho[l] = 0.125;
	  by[l]  = -1.0;
	}
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif

}
