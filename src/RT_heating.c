//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Heating_cpu () {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Vx);
  INPUT(Vy);
  INPUT(Vz);
  OUTPUT(QPlus);
#ifdef IRRADIATION
  INPUT (StellarRad);
#endif
//<\USER_DEFINED>


//<EXTERNAL>
  real* rho       = Density->field_cpu;
  real* vx        = Vx->field_cpu;
  real* vy        = Vy->field_cpu;
  real* vz        = Vz->field_cpu;
  real* qplus     = QPlus->field_cpu;
#ifdef IRRADIATION
  real* stellarrad = StellarRad->field_cpu;
#endif
  int pitch       = Pitch_cpu;
  int stride      = Stride_cpu;
  int size_x      = XIP;
  int size_y      = Ny+2*NGHY-1;
  int size_z      = Nz+2*NGHZ-1;
  real dx         = Dx;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real div_v;
  real vror;
  real heat;
  real cotr;
//<\INTERNAL>
 
//<CONSTANT>
// real NU(1);
// real Sxj(Ny+2*NGHY);
// real Syj(Ny+2*NGHY);
// real Szj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syk(Nz+2*NGHZ);
// real Szk(Nz+2*NGHZ);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real InvVj(Ny+2*NGHY);
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
//Evaluate centered divergence.
	div_v = 0.0;
#ifdef X
	div_v += (vx[lxp]-vx[l])*SurfX(j,k);
#endif
#ifdef Y
	div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
#endif
#ifdef Z
	div_v += (vz[lzp]*SurfZ(j,k+1)-vz[l]*SurfZ(j,k));
#endif
	div_v *= 2.0/3.0*InvVol(j,k);

	heat  = 2.0*pow((vy[lyp]-vy[ll])/zone_size_y(j,k),2.0);
	vror  = (vy[ll]+vy[lyp])/(ymin(j)+ymin(j+1));
	cotr  = cos(zmed(k))/sin(zmed(k))/ymed(j);
	heat += 2.0*pow((vz[lzp]-vz[ll])/zone_size_z(j,k)+vror,2.0);
	heat += 2.0*pow((vx[lxp]-vx[ll])/zone_size_x(j,k)+vror+	\
			.5*(vz[ll]+vz[lzp])*cotr,2.0);
	heat -= 2.0/3.0*div_v*div_v;
	heat += pow(.25*((vz[lyp]+vz[lzp+pitch])-(vz[lym]+vz[lzp-pitch]))/zone_size_y(j,k)- \
		    (vz[ll]+vz[lzp])/(ymed(j+1)+ymed(j))+\
		    .25*((vy[lzp]+vy[lzp+pitch])-(vy[lzm]+vy[lzm+pitch]))/zone_size_z(j,k),2.0);
	heat += pow(.25*((vx[lzp]+vx[lxp+stride])-(vx[lzm]+vx[lxp-stride]))/zone_size_z(j,k)- \
		    cotr*.5*(vx[ll]+vx[lxp])+				\
		    .25*((vz[lxp]+vz[lxp+stride])-(vz[lxm]+vz[lxm+stride]))/zone_size_x(j,k),2.0);
	heat += pow(.25*((vy[lxp]+vy[lxp+pitch])-(vy[lxm]+vy[lxm+stride]))/zone_size_x(j,k)+ \
		    .25*((vx[ll+pitch]+vx[lxp+pitch])-(vx[ll-pitch]+vx[lxp-pitch]))/zone_size_y(j,k)- \
		    .5*(vx[ll]+vx[lxp])/ymed(j),2.0);
	qplus[ll] = NU*heat*rho[ll];
#ifdef IRRADIATION
	qplus[ll] += stellarrad[ll];
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
