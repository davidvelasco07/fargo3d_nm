#include "fargo3d.h"
void boundaries() {
  //Current_Fluid is the global fluid rank
  //It is asigned in SelectFluid
  if (!PERIODICZ) {
#ifdef Z
    if(LocalBC[_Z_][0] > 0)
      boundary_zmin[Current_Fluid]();
    if(LocalBC[_Z_][1] > 0)
      boundary_zmax[Current_Fluid]();
#endif
  }
  if (!PERIODICY) {
#ifdef Y
    if(LocalBC[_Y_][0] > 0)
      boundary_ymin[Current_Fluid]();
    if(LocalBC[_Y_][1] > 0)
      boundary_ymax[Current_Fluid]();
#endif
  }
  if (!PERIODICX) {
#ifdef X
    if(LocalBC[_X_][0] > 0)
      boundary_xmin[Current_Fluid]();
    if(LocalBC[_X_][1] > 0)
      boundary_xmax[Current_Fluid]();
#endif
  }
}

