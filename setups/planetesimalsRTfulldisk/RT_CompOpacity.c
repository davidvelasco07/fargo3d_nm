//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ComputeOpacity_cpu () {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Temperature);
  OUTPUT(OpaPlanck);
//<\USER_DEFINED>


//<EXTERNAL>
  real* opa = OpaPlanck->field_cpu;
  real* t   = Temperature->field_cpu;
  real* dens= Density->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real kappa = KAPPA;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real temp;
  real rho;
  real power1;
  real power2;
  real power3;
  real t234;
  real t456;
  real t678;
  real ak1;
  real ak2;
  real ak3;
  real bk3;
  real bk4;
  real bk5;
  real bk6;
  real bk7;
  real ts4;
  real rho13;
  real rho23;
  real ts42;
  real ts44;
  real ts48;
  real o3;
  real o4;
  real o5;
  real o4an;
  real o3an;
  real t2;
  real t4;
  real t8;
  real t10;
  real o1;
  real o2;
  real o1an;
  real o2an;
  real o6;
  real o7;
  real o6an;
  real o7an;
  real opacity;
//<\INTERNAL>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
//	ts4 = 0.0;
//	ts42 = 0.0;
//	opacity = 0.0;
	ll = l;
//
//	// In the module below we need the temperature and density in cgs units
//	temp = t[ll]/(G*MSTAR/(R0*R_MU))*(G_CGS*MSTAR_CGS/(R0_CGS*R_MU_CGS));
//	rho  = dens[ll]/(MSTAR/(R0*R0*R0))*(MSTAR_CGS/(R0_CGS*R0_CGS*R0_CGS));
//
//	//  The subroutine below has been taken from the JUPITER code,
//	//  itself inherited from FARGOCA, probably itself inherited
//	//  of another code (FM, 07/06/2014). The units MUST be cgs
//	//  for this part to work, in its present state.
//
//	// Opacity routine: Doug Lin's Opacity Function (modified)
//
//	//   data power1,power2,power3/4.44444444e-2, 2.381e-2, 2.267e-1/
//	power1=4.44444444e-2;
//	power2=2.381e-2;
//	power3=2.267e-1;
//
//	//    data t234,t456,t678/1.6e3, 5.7e3, 2.28e6/
//	t234=1.6e3;
//	t456=5.7e3;
//	t678=2.28e6;
//
//	//  coefficients for opacity laws 1, 2, and 3 in cgs units.
//	//  data ak1,ak2,ak3/2.e-4, 2.e16, 5.e-3/
//	ak1=2.e-4;
//	ak2=2.e16;
//	ak3=5.e-3;
//
//	/*   coefficients for opacity laws 3, 4, 5, 6, 7, and 8 in T_4 units.
//	     data bk3,bk4,bk5,bk6,bk7   50., 2.e-2, 2.e4, 1.d4, 1.5d10 */
//	bk3 = 50.;
//	bk4 = 2.e-2;
//	bk5 = 2.e4;
//	bk6 = 1.e4;
//	bk7 = 1.5e10;
//
//	/*  test T against (T_23 * T_34 * T_34)**0.333333333 */
//
//	if( temp < t234*pow(rho,power1)){
//	  //   different powers of temperature
//	  t2=temp*temp;
//	  t4=t2*t2;
//	  t8=t4*t4;
//	  t10=t8*t2;
//	  //   disjoint opacity laws
//	  o1=ak1*t2;
//	  o2=ak2*temp/t8;
//	  o3=ak3*temp;
//	  //   parameters used for smoothing
//	  o1an=o1*o1;
//	  o2an=o2*o2;
//	  //   smoothed and continuous opacity law for regions 1, 2, and 3.
//	  opacity=pow(pow(o1an*o2an/(o1an+o2an),2.0)+pow(o3/(1.0+1.e22/t10),4.0),0.25);
//	}
//  
//	if (temp > t234*pow(rho,power1) && temp < t456*pow(rho,power2)){
//	  
//	  //   to avoid overflow
//	  ts4=1.e-4*temp;
//	  rho13=pow(rho,0.333333333);
//	  rho23=rho13*rho13;
//	  ts42=ts4*ts4;
//	  ts44=ts42*ts42;
//	  ts48=ts44*ts44;
//	  
//	  //   disjoint opacity laws for 3, 4, and 5.
//	  o3=bk3*ts4;
//	  o4=bk4*rho23/(ts48*ts4);
//	  o5=bk5*rho23*ts42*ts4;
//	  //   parameters used for smoothing
//	  o4an=pow(o4,4.0);
//	  o3an=pow(o3,4.0);
//	  //  smoothed and continuous opacity law for regions 3, 4, and 5.
//	  opacity=pow((o4an*o3an/(o4an+o3an))+pow(o5/(1.0+6.561e-5/ts48),4.0),0.25);
//	}
//        
//	if (temp > t456*pow(rho,power2)){
//	  if (temp < t678*pow(rho,power3) || rho <= 1.0e-10){
//	    //  to avoid overflow
//	    ts4=1.e-4*temp;
//	    rho13=pow(rho,0.333333333);
//	    rho23=rho13*rho13;
//	    ts42=ts4*ts4;
//	    ts44=ts42*ts42;
//	    ts48=ts44*ts44;
//	    
//	    //   disjoint opacity laws for 5, 6, and 7.
//	    o5=bk5*rho23*ts42*ts4;
//	    o6=bk6*rho13*ts44*ts44*ts42;
//	    o7=bk7*rho/(ts42*sqrt(ts4));
//	    //   parameters used for smoothing
//	    o6an=o6*o6;
//	    o7an=o7*o7;
//	    //   smoothed and continuous opacity law for regions 5, 6, and 7.
//	    opacity=pow(pow(o6an*o7an/(o6an+o7an),2.0)+\
//			pow(o5/(1.0+pow(ts4/(1.1*pow(rho,0.04762)),10.0)),4.0),0.25);
//	  }
//	  else {
//	    //  no scattering!
//	    opacity=bk7*rho/(ts42*sqrt(ts4));
//	  }
//	}
//	// The opacity thus obtained was in cgs units (cm^2/g). We
//	// convert it back to the current units
	opacity = kappa;
	opa[ll] = opacity/(R0_CGS*R0_CGS/MSTAR_CGS)*R0*R0/MSTAR;
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
