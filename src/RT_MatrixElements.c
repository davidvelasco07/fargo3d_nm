//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ComputeMatrixElements_cpu (real dt) {

//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Energyrad);
  INPUT(Density);
  INPUT(DiffCoef);
  INPUT(Temperature);
  INPUT(OpaPlanck);
  INPUT(Eta1);
  INPUT(Eta2);
  OUTPUT(LeftX);
  OUTPUT(LeftY);
  OUTPUT(LeftZ);
  OUTPUT(RightX);
  OUTPUT(RightY);
  OUTPUT(RightZ);
  OUTPUT(RHS);
  OUTPUT(BB);
//<\USER_DEFINED>


//<EXTERNAL>
  real* ener     = Energy->field_cpu;
  real* enrad    = Energyrad->field_cpu;
  real* dens     = Density->field_cpu;
  real* diffcoef = DiffCoef->field_cpu;
  real* temper   = Temperature->field_cpu;
  real* opap     = OpaPlanck->field_cpu;
  real* eta1     = Eta1->field_cpu;
  real* eta2     = Eta2->field_cpu;
  real* leftx    = LeftX->field_cpu;
  real* lefty    = LeftY->field_cpu;
  real* leftz    = LeftZ->field_cpu;
  real* rightx   = RightX->field_cpu;
  real* righty   = RightY->field_cpu;
  real* rightz   = RightZ->field_cpu;
  real* rhs      = RHS->field_cpu;
  real* bb       = BB->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP;
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
  real dx        = Dx;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real temper3;
//<\INTERNAL>
  
//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>
  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=1; k<size_z; k++) {
#endif
#ifdef Y
    for (j=1; j<size_y; j++) {
#endif
#ifdef X
      for (i=XIM; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
	rightx[ll]  = -.5*dt*(diffcoef[ll]+diffcoef[lxp])/(zone_size_x(j,k)*zone_size_x(j,k));
	leftx[ll]   = -.5*dt*(diffcoef[ll]+diffcoef[lxm])/(zone_size_x(j,k)*zone_size_x(j,k));
	righty[ll]  = -.5*dt*(diffcoef[ll]+diffcoef[lyp])*\
	  pow(ymin(j+1)/(ymed(j)*(ymin(j+1)-ymin(j))),2.0);
       	lefty[ll]   = -.5*dt*(diffcoef[ll]+diffcoef[lym])*\
	  pow(ymin(j)/(ymed(j)*(ymin(j+1)-ymin(j))),2.0);
	rightz[ll]  = -.5*dt*(diffcoef[ll]+diffcoef[lzp])*\
	  sin(zmin(k+1))/((cos(zmin(k))-cos(zmin(k+1)))*(ymed(j)*zone_size_z(j,k)));
        leftz[ll]   = -.5*dt*(diffcoef[ll]+diffcoef[lzm])*\
	  sin(zmin(k))/((cos(zmin(k))-cos(zmin(k+1)))*(ymed(j)*zone_size_z(j,k)));
	bb[ll]      = -leftx[ll]-rightx[ll]-lefty[ll]-righty[ll]-leftz[ll]-rightz[ll];
	temper3     = pow(temper[ll],3.0);
	bb[ll]      = 1.0+dt*dens[ll]*opap[ll]*(C0-16.0*STEFANK*temper3*eta2[ll])+bb[ll];
	rhs[ll]     = enrad[ll]+4.0*dt*dens[ll]*opap[ll]*STEFANK*temper3*(4.0*eta1[ll]-3.0*temper[ll]);
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
