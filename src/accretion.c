//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void compute_accretion(real dt) {
  int i,j,k,n;
  real* xp = Sys->x;
  real* yp = Sys->y;
  real* zp = Sys->z;
  real* vxp = Sys->vx;
  real* vyp = Sys->vy;
  real* vzp = Sys->vz;
  real* mp = Sys->mass;
  int nb = Sys->nb;
  real rhop = 3*(MSTAR/(R0*R0*R0))/(MSTAR_CGS/(R0_CGS*R0_CGS*R0_CGS)); //3g/cm^3 -> c.u
  real Rp,dM;

  for (n=0;n<nb;n++){
    Accretion(xp[n], yp[n], zp[n], vxp[n], vyp[n], vzp[n], mp[n], dt);
    dM = reduction_full_SUM(Energy, NGHY, Ny+NGHY, 0, Nz+2*NGHZ);
    M_acc[FluidIndex][n] += dM;
    Rp = pow(3*mp[n]/(4*M_PI*rhop),1./3);
    Lum[FluidIndex][n] = G*mp[n]*dM/dt/Rp; 
  }
}

void Accretion_cpu (real xp, real yp, real zp, real vxp, real vyp, real vzp, real Mp, real dt) {
//<USER_DEFINED>
  INPUT(Density);
  OUTPUT(Density);
  OUTPUT(Energy);
#ifdef X
  INPUT(Vx);
#endif
#ifdef Y
  INPUT(Vy);
#endif
#ifdef Z
  INPUT(Vz);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho_d = Density->field_cpu;
#ifdef X 
  real* vx_d = Vx->field_cpu;
#endif
#ifdef Y
  real* vy_d = Vy->field_cpu;
#endif
#ifdef Z
  real* vz_d = Vz->field_cpu;
#endif
  real* accreted = Energy->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real stokes = 1./Coeffval[0];
//<\EXTERNAL>
  
//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int ll;
  int ip;
  int jp;
  int kp;
  int lp;
  real x;
  real y;
  real z;
  real vx;
  real vy;
  real vz;
  real dx;
  real dy;
  real dz;
  real dxp;
  real dyp;
  real dzp;
  real r_polar;
  real dist;
  real frac;
  real removed;
  real v_hw;
  real Mdot;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real Sxj(Ny+2*NGHY);
// real Syj(Ny+2*NGHY);
// real Szj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syk(Nz+2*NGHZ);
// real Szk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
// real OMEGAFRAME(1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;
#ifdef Z
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=1; j<size_y; j++) {
#endif
#ifdef X
      for(i=1; i<size_x; i++) {
#endif
//<#>
	ll = l;
  r_polar = sqrt(xp*xp+yp*yp);
#ifdef CARTESIAN
  x = xp;
  y = yp;
  z = zp;
#endif
#ifdef CYLINDRICAL
  x = atan2(yp,xp);
  y = sqrt(xp*xp+yp*yp);
  z = zp;
#endif
#ifdef SPHERICAL
  x = atan2(yp,xp);
  y = sqrt(xp*xp+yp*yp+zp*zp);
  z = atan2(sqrt(xp*xp+yp*yp),zp);
  //v_phi
  vx = (xp*vyp-vxp*yp)/r_polar;
  //v_r
  vy = (xp*vxp + yp*vyp + zp*vzp)/y;
  //v_theta
  vz = ((xp*vxp + yp*vyp)*zp/r_polar - r_polar*vzp)/y;
#endif
  dx = xmed(1)-xmed(0);
  dy = ymed(1)-ymed(0);
  dz = zmed(1)-zmed(0);
  
  //Indices of cell hosting the planet
  //ip = floor((x-xmed(0))/dx+.5);
  //jp = floor((y-ymed(0))/dy+.5);
  //kp = floor((z-zmed(0))/dz+.5);
  //lp = ip + jp*pitch + kp*stride;
  //This should be an average over the nearest cells to the planet
  //v_hw = sqrt(pow((vx_d[lp]+OMEGAFRAME*ymed(j)*sin(zmed(k)))-vx,2)+pow(vy_d[lp]-vy,2)+pow(vz_d[lp]-vz,2));
  //RB = G*Mp/v_hw/v_hw;

  dxp=fabs(x-xmed(i));
	dyp=fabs(y-ymed(j));
	dzp=fabs(z-zmed(k));
  //drp=sqrt(dxp*dxp+dyp*dyp+dzp*dzp);
  
  frac=1;
  //Just in case let's clean whatever was in this field
  accreted[ll] = 0;
	if(dxp<dx && dyp<dy && dzp<dz){
      //Mdot = rho_d[ll]*M_PI*RB*RB/v_hw;
	    frac *= 1-dxp/dx;
	    frac *= 1-dyp/dy;
	    frac *= 1-dzp/dz;
      //removed = frac*Mdot*dt;
      removed = frac*rho_d[ll];
      rho_d[ll] -= removed;
      accreted[ll] = Vol(j,k)*removed;
      //printf("(%d,%d,%d)->%g, dens=%lf, removed=%lf, Vol=%.18g, Macc=%.18g\n",i,j,k,frac,rho_d[ll],removed,Vol(j,k),accreted[ll]);
  }
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
