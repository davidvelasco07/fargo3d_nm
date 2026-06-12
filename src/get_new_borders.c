#include <fargo3d.h>

void GetNewBorders() {
  StateVector v;
  OrbitalElements o;

  real x,y,r,a;

  int i;

  real amin=1.0e+30*R0;
  real amax=0.0;
  
  for (i=0; i < Sys->nb; i++) {

    v.x = Sys->x[i];
    v.y = Sys->y[i];
    v.z = Sys->z[i];
    v.vx = Sys->vx[i];
    v.vy = Sys->vy[i];
    v.vz = Sys->vz[i];
    o = SV2OE (v,MSTAR+Sys->mass[i]);

    if (amax < o.a) {
      amax = o.a;
    }
    if (amin > o.a) {
      amin = o.a;
    }

  }
  
  InnerBorder = amin * pow(RESONANCE,-2.0/3.0);
  OuterBorder = amax * pow(RESONANCE, 2.0/3.0);
  
}
