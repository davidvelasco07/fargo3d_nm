#include <fargo3d.h>

void RadialScan (Field *array) { //Integrates array (\tau, generally) radially (radial prefix sum)

  FARGO_SAFE (InProcessRadialScan (array)); // cpu or gpu : radial scan
					    // of a CPU block, radial
					    // ghosts excluded
  FARGO_SAFE (GetBackgroundSlice (array)); // standard C2CUDA
  FARGO_SAFE (InterCPUScan ()); //not standard C2CUDA
  FARGO_SAFE (OffsetPatches (array)); // cpu or gpu; standard C2CUDA
  FARGO_SAFE (SetBackgroundSlice (array)); // standard C2CUDA
}
