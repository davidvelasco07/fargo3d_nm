#include "fargo3d.h"

real Resistivity (real y, real z) {
  //The user can freely redefine this function
  //By defining his own function, in a file named
  //resistivity.c, in his own setup directory,
  //which will be built instead of this file.

  return NU/PRANDTL;
}

