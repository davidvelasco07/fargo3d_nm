#include "fargo3d.h"

void _CondInit() {
  
  int i,j,k;
  real r, omega;
  
  real *rho  = Density->field_cpu;
  real *cs   = Energy->field_cpu;
  real *vphi = Vx->field_cpu;
  real *vr   = Vy->field_cpu;
  
  real rhog, rhod;
  real vk;
  
  i = j = k = 0;
  
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {
	
	r     = Ymed(j);
	omega = sqrt(G*MSTAR/r/r/r);                       //Keplerian frequency
	rhog  = SIGMA0*pow(r/R0,-SIGMASLOPE);              //Gas surface density
        rhod  = rhog*EPSILON;                              //Dust surface density

	if (Fluidtype == GAS) {
	  rho[l]   = rhog;
	  vphi[l]  = omega*r*sqrt(1.0 + pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
				  (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));
	  vr[l]    = 0.0;
	  cs[l]    = ASPECTRATIO*pow(r/R0,FLARINGINDEX)*omega*r;
	}
	
	if (Fluidtype == DUST) {
	  rho[l]  = rhod;
	  vphi[l] = omega*r;
	  vr[l]   = 0.0;
	  cs[l]   = 0.0;
	}
	
	vphi[l] -= OMEGAFRAME*r;
	
      }
    }
  }
}

void CondInit() {

  int feedback = YES;
  char dust_name[MAXNAMELENGTH];
  int global_index;
  int id;
  
  printf("fluid index = %d \n", FluidColor);
  
  for (id = 0; id<NFluids_per_rank; id++) {
    global_index = FluidColor*NFluids_per_rank+id;
    if (global_index == 0) {
      Fluids[id] = CreateFluid("gas",GAS);
      SelectFluid(id);
      _CondInit();
    }
    else {
      sprintf(dust_name,"dust%d",global_index);
      Fluids[id]  = CreateFluid(dust_name, DUST);
      SelectFluid(id);
      _CondInit();
    }
  }


}
