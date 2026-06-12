#include "fargo3d.h"
void boundaries() {
  if (!PERIODICZ) {
#ifdef Z
    if (LocalBC[_Z_][0] > 0)
      FARGO_SAFE(boundary_zmin());
    if (LocalBC[_Z_][1] > 0)
      FARGO_SAFE(boundary_zmax());
#endif
  }
  if (!PERIODICY) {
#ifdef Y
    if (LocalBC[_Y_][0] > 0){
      FARGO_SAFE(boundary_ymin());
    }
    if (LocalBC[_Y_][1] > 0)
      FARGO_SAFE(boundary_ymax());
#endif
  }
  if (!PERIODICX) {
#ifdef X
    if (LocalBC[_X_][0] > 0)
      FARGO_SAFE(boundary_xmin());
    if (LocalBC[_X_][1] > 0)
      FARGO_SAFE(boundary_xmax());
#endif
  }
}

