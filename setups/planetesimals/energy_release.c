#include "fargo3d.h"

void energy_release (real dt) {
  real lum, mpl, xplanet, yplanet, zplanet;
  real azim, rad;
  int ip, jp, kp;
  xplanet = Sys->x_cpu[0];
  yplanet = Sys->y_cpu[0];
  zplanet = Sys->z_cpu[0];
  mpl     = Sys->mass_cpu[0];
  lum = G*pow(mpl,5./3.)*pow(4.*M_PI/3.*RHOSOLID,1./3.)/DOUBLINGTIME;
  azim = atan2 (yplanet, xplanet);
  rad  = sqrt (xplanet*xplanet+yplanet*yplanet+zplanet*zplanet);
  ip = (int)((azim-XMIN)/(XMAX-XMIN)*(real)NX);
  jp = (int)((rad-YMIN)/(YMAX-YMIN)*(real)NY);
  kp = 0;
  energy_release_b (lum, dt, ip, jp, kp, y0cell, z0cell);
}
