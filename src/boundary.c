#include "fargo3d.h"
void boundaries() {

  int global_index = FluidColor*NFluids_per_rank+FluidIndex;
  
  if (!PERIODICZ) {
#ifdef Z
    if(LocalBC[_Z_][0] > 0)
      boundary_zmin[global_index]();
    if(LocalBC[_Z_][1] > 0)
      boundary_zmax[global_index]();
#endif
  }
  if (!PERIODICY) {
#ifdef Y
    if(LocalBC[_Y_][0] > 0)
      boundary_ymin[global_index]();
    if(LocalBC[_Y_][1] > 0)
      boundary_ymax[global_index]();
#endif
  }
  if (!PERIODICX) {
#ifdef X
    if(LocalBC[_X_][0] > 0)
      boundary_xmin[global_index]();
    if(LocalBC[_X_][1] > 0)
      boundary_xmax[global_index]();
#endif
  }
}

