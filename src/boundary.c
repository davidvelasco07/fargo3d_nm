#include "fargo3d.h"
void boundaries() {

  int global_index = FluidColor*NFluids_per_rank+FluidIndex;
  
  if (!PERIODICZ) {
#ifdef Z
    if(Gridd.bc_down)
      boundary_zmin[global_index]();
    if(Gridd.bc_up)
      boundary_zmax[global_index]();
#endif
  }
  if (!PERIODICY) {
#ifdef Y
    if(Gridd.bc_left)
      boundary_ymin[global_index]();
    if(Gridd.bc_right)
      boundary_ymax[global_index]();
#endif
  }
#ifdef GHOSTSX 
  Fill_GhostsX();
#endif
}

