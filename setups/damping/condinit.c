#include "fargo3d.h"

void _CondInit(real rho0, real v1, real v2, real v3){
  
  int i,j,k;
  
  real* rho = Density->field_cpu;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
  real* vz = Vz->field_cpu;
  real* e  = Energy->field_cpu;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j = 0; j<Ny+2*NGHY; j++) {
      for (i = 0; i<Nx; i++) {
	rho[l] = rho0;
	vx[l] = v1;
	vy[l] = v2;
	vz[l] = v3;
	e[l]  = 0.0;
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
      _CondInit(1.0e3, 0.0,2.0,0.0);
    }
    else {
      sprintf(dust_name,"dust%d",global_index);
      Fluids[id]  = CreateFluid(dust_name, DUST);
      SelectFluid(id);
      _CondInit(1.0e-3,0.0,-3.0*global_index,0.0);
    }
  }

}

