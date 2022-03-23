//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void compute_accretion(real dt) {
  int i,j,k,n;
  real* xp = Sys->x_cpu;
  real* yp = Sys->y_cpu;
  real* zp = Sys->z_cpu;
  real* vxp = Sys->vx;
  real* vyp = Sys->vy;
  real* vzp = Sys->vz;
  real* mp = Sys->mass;
  real M_acc[10];
  char filename[MAXLINELENGTH];
  FILE *Out;
  int nb = Sys->nb;

  for (n=0;n<nb;n++){
    Accretion(xp[n], yp[n], zp[n], vxp[n], vyp[n], vzp[n], mp[n], dt);
    M_acc[n] = reduction_full_SUM(Energy, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
  }
  if (CPU_Master)
  {
    sprintf(filename, "%smonitor/%s/M_acc.dat", OUTPUTDIR, Fluids[FluidIndex]->name);
    Out = fopen_prs(filename, "a+");
    fprintf(Out, "%.12g\t%.12g\n", PhysicalTime, M_acc[0]);
    fclose(Out);
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
  real r_polar;
  real dist;
  real v_hw;
  real M_Ormel;
  real b;
  real t_enc;
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
	dx = xmed(i)-xp;
	dy = ymed(j)-yp;
	dz = zmed(k)-zp;
  x = xp;
  y = yp;
  z = zp;
  vx = vxp;
  vy = vyp;
  vz = vzp;
#endif
#ifdef CYLINDRICAL
	dx = ymed(j)*cos(xmed(i))-xp;
	dy = ymed(j)*sin(xmed(i))-yp;
	dz = zmed(k)-z;
  x = atan2(yp,xp);
  y = sqrt(xp*xp+yp*yp);
  z = zp;
  //v_phi
  vx = (xp*vyp-vxp*yp)/y + OMEGAFRAME*ymed(j);
  //v_r
  vy = (xp*vxp+yp*vyp)/y;
  //v_z
  vz = vzp;
#endif
#ifdef SPHERICAL
	dx = ymed(j)*cos(xmed(i))*sin(zmed(k))-xp;
	dy = ymed(j)*sin(xmed(i))*sin(zmed(k))-yp;
	dz = ymed(j)*cos(zmed(k))-zp;
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
  //Indices of cell hosting the planet
  ip = (x-xmed(0))/(xmed(1)-xmed(0));
  jp = (y-ymed(0))/(ymed(1)-ymed(0));
  kp = (z-zmed(0))/(zmed(1)-zmed(0));
  lp = ip + jp*pitch + kp*stride;
#ifdef ACCRETION
  //Headwind velocity
  v_hw = sqrt(pow(vx_d[lp]+OMEGAFRAME*ymed(jp)*sin(zmed(kp))-vx,2)+pow(vy_d[lp]-vy,2)+pow(vz_d[lp]-vz,2));
  //printf("vp=(%lf,%lf,%lf), vd=(%lf,%lf,%lf), vhw=%lf\n",vx,vy,vz,vx_d[lp]+OMEGAFRAME*ymed(jp)*sin(zmed(kp)),vy_d[lp],vz_d[lp],v_hw);
  //printf("v_phi=%lf, OmegaFrame=%lf,r=%lf, sin(theta)=%lf\n",vx_d[lp],OMEGAFRAME,ymed(jp),sin(zmed(kp)));
  M_Ormel = pow(v_hw,3)/(8*G*OMEGAFRAME*stokes)/MSTAR;
  if(Mp > M_Ormel){
    //Shear regime
    b = pow(G*Mp*stokes,1./3);
    t_enc = 1/OMEGAFRAME;
  }
  else{
    //Headwind regime
    b = sqrt(2*G*Mp*OMEGAFRAME*stokes/v_hw);
    t_enc = 2*b/v_hw;
  }
  b = MIN(b, y*pow(Mp/MSTAR/3.0,1./3.));
  b = y*pow(Mp/MSTAR/3.0,1./3.);
  //printf("v_hw=%lf, M_Ormel=%lf, b=%lf\n",v_hw, M_Ormel, b);
  t_enc = 1.01*dt;//MAX(1.01*dt,t_enc);
  dist = sqrt(dx*dx+dy*dy+dz*dz);
  
  //if(dist < b){
  if ((i == ip) && (j == jp) && (k == kp)) {
    accreted[ll] = dt/t_enc*rho_d[ll];
    rho_d[ll] -= accreted[ll];
    accreted[ll] = Vol(j,k)*accreted[ll];//*(1-hidden[ll])
  }
#endif
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
