//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

//#define BLOCK_X 8
//#define BLOCK_Y 8
//#define BLOCK_Z 1

// Fills forces arrays of the force experienced by the planet

void _ComputeForce_cpu(real x, real y, real z, real rsmoothing, real mass) {

//<USER_DEFINED>
  INPUT(Density);
  OUTPUT(Mmx);
  OUTPUT(Mpx);
  OUTPUT(Mmy);
  OUTPUT(Mpy);
  OUTPUT(Vx_temp);
  OUTPUT(Vy_temp);
//<\USER_DEFINED>

//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* fxi  = Mmx->field_cpu;
  real* fyi  = Mpx->field_cpu;
  real* fzi  = Mmy->field_cpu;
  real* fxhi = Mpy->field_cpu;
  real* fyhi = Vx_temp->field_cpu;
  real* fzhi = Vy_temp->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+NGHY;
  int size_z = Nz+NGHZ;
  real rh = pow(mass/3./MSTAR, 1./3.)*sqrt(x*x+y*y);
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  real cellmass;
  real distance;
  real dist2;
  real dx;
  real dy;
  real dz;
  real InvDist3; 
  real hill_cut; 
  real planet_distance;
//<\INTERNAL>

//<CONSTANT>
// real Sxj(Ny+2*NGHY);
// real Syj(Ny+2*NGHY);
// real Szj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syk(Nz+2*NGHZ);
// real Szk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=NGHZ; k<size_z; k++) {
#endif
#ifdef Y
    for (j = NGHY; j<size_y; j++) {
#endif
#ifdef X
      for (i = 0; i<size_x; i++) {
#endif

//<#>
	/* Be careful! Vol(j,k) must be the right 
	   expresion for 2D simulations also (ie: surface area). */

	cellmass = Vol(j,k)*dens[l];
#ifdef CYLINDRICAL
	dx = ymed(j)*cos(xmed(i))-x;
	dy = ymed(j)*sin(xmed(i))-y;
	dz = 0.0;
#endif
	dist2 = dx*dx+dy*dy;
	/* New default exclusion function */

#ifdef HILLCUT 
	planet_distance=sqrt(dist2);

	if (planet_distance/rh < 0.5)
	  hill_cut = 0.0;
	else {
	  if (planet_distance > rh)
	    hill_cut = 1.0;
	  else
	    hill_cut = pow(sin((planet_distance/rh-.5)*M_PI),2.);
	}
#endif	

	dist2 += rsmoothing*rsmoothing;
	distance = sqrt(dist2);
	InvDist3 = 1.0/(dist2*distance);
	InvDist3 *= G*cellmass;
	
	fxi[l]  = dx*InvDist3;
	fyi[l]  = dy*InvDist3;
	fzi[l]  = 0.0;
#ifdef HILLCUT
	fxhi[l] = dx*InvDist3*hill_cut;
	fyhi[l] = dy*InvDist3*hill_cut;
	fzhi[l] = 0.0;
#else
	fxhi[l] = 0.0;
	fyhi[l] = 0.0;
	fzhi[l] = 0.0;
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

//<LAST_BLOCK>
  
  int index;
  real aa = sqrt(x*x+y*y+z*z);

  /*This part must be improved for 3D*/
  index = (int)(((aa-Ymin(NGHY))/(Ymin(Ny+NGHY+1)-Ymin(NGHY)))*(real)Ny) + NGHY;
  
  if (index >= NGHY) {
    if(index < Ny+NGHY) {
      /*Inner Force*/
      localforce[0]  = reduction_full_SUM(Mmx, NGHY, index, NGHZ, Nz+NGHZ);
      localforce[1]  = reduction_full_SUM(Mpx, NGHY, index, NGHZ, Nz+NGHZ);
      localforce[2]  = reduction_full_SUM(Mmy, NGHY, index, NGHZ, Nz+NGHZ);
      localforce[3]  = reduction_full_SUM(Mpy, NGHY, index, NGHZ, Nz+NGHZ);
      localforce[4]  = reduction_full_SUM(Vx_temp, NGHY, index, NGHZ, Nz+NGHZ);
      localforce[5]  = reduction_full_SUM(Vy_temp, NGHY, index, NGHZ, Nz+NGHZ);
      /*Outer Force*/
      localforce[6]  = reduction_full_SUM(Mmx, index, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[7]  = reduction_full_SUM(Mpx, index, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[8]  = reduction_full_SUM(Mmy, index, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[9]  = reduction_full_SUM(Mpy, index, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[10] = reduction_full_SUM(Vx_temp, index, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[11] = reduction_full_SUM(Vy_temp, index, Ny+NGHY, NGHZ, Nz+NGHZ);
    }
    /*All is Inner Force*/
    else{
      localforce[0]   = reduction_full_SUM(Mmx, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[1]   = reduction_full_SUM(Mpx, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[2]   = reduction_full_SUM(Mmy, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[3]   = reduction_full_SUM(Mpy, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[4]   = reduction_full_SUM(Vx_temp, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[5]   = reduction_full_SUM(Vy_temp, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
      localforce[6]   = 0.0;
      localforce[7]   = 0.0;
      localforce[8]   = 0.0;
      localforce[9]   = 0.0;
      localforce[10]  = 0.0;
      localforce[11]  = 0.0;
    }
  }
  /*All is Outer Force*/
  else{
    localforce[0]  = 0.0;
    localforce[1]  = 0.0;
    localforce[2]  = 0.0;
    localforce[3]  = 0.0;
    localforce[4]  = 0.0;
    localforce[5]  = 0.0;
    localforce[6]  = reduction_full_SUM(Mmx, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
    localforce[7]  = reduction_full_SUM(Mpx, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
    localforce[8]  = reduction_full_SUM(Mmy, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
    localforce[9]  = reduction_full_SUM(Mpy, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
    localforce[10] = reduction_full_SUM(Vx_temp, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
    localforce[11] = reduction_full_SUM(Vy_temp, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
  }

//<\LAST_BLOCK>

}
