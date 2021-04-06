#include "fargo3d.h"

void _CondInit(int id) {
  int i,j,k, n;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *rho;
  
  real omega, rho_s, alp;
  real r, r3, xi, beta, h;
  
  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;

  real stokes_plus[NFLUIDS];
  real stokes[NFLUIDS-1];
  real epsilons[NFLUIDS-1];
  real sq = SQ;
  real slope, hd, sigma;
  // Stokes numbers                                                                       
  real smax = SMAX;
  real smin = SMIN;
  real ds   = (log(smax)-log(smin))/(NFLUIDS-1);
  for(n=0;n<NFLUIDS;n++){
    stokes_plus[n] = smin*exp(ds*n);
  }

  //Density distribution                                                                  
  slope = 4-sq;
  for(n=0; n<NFLUIDS-1; n++){
    if(slope != 0.) {
      epsilons[n]  = pow(stokes_plus[n+1],slope) - pow(stokes_plus[n],slope) ;
      epsilons[n] *= EPSILON/(pow(smax, slope) - pow(smin,slope));
    }
    else{
      epsilons[n]  = log(stokes_plus[n+1]/stokes_plus[n]);
      epsilons[n] *= EPSILON/log(smax/smin);
    }
    stokes[n] = stokes_plus[n+1];
    if( NFLUIDS == 2) stokes[n] = SMAX;
  }
  
  sigma = SIGMA0;
  if(id > 0) {
    sigma   = SIGMA0*epsilons[id-1];
#ifdef STOKESNUMBER
    Coeffval[0]   = 1.0/stokes[id-1];
#endif
#ifdef DUSTSIZE
    Coeffval[1]   = 1.0/(stokes[id-1]*R0/R0_CGS);    
    Coeffval[2]   = RHOSOLID/(MSTAR_CGS/(R0_CGS*R0_CGS*R0_CGS))*(MSTAR/(R0*R0*R0));
#endif
    //if(CPU_Master) printf("%1.16f \t %1.16f \n", stokes[id-1], epsilons[id-1]);
  }


  if(Fluidtype == DUST) hd = MIN( 1.0, sqrt( (DELTA+NU)/( stokes[id-1] + (DELTA+NU)) ));
  if(Fluidtype == GAS)  hd = 1.0;
  
  xi = SIGMASLOPE+1.+FLARINGINDEX;
  beta = 1.-2*FLARINGINDEX;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      r = Ymed(j);
      r3 = r*r*r;
      omega = sqrt(G*MSTAR/(r3));

      h = ASPECTRATIO*pow(r/R0,FLARINGINDEX);

      for (i=NGHX; i<Nx+NGHX; i++) {
	//Dust enery fills with gas energy for the dragforce coeffcient calculation
	e[l] = h*sqrt(G*MSTAR/r);

	v2[l] = v3[l] = 0.0;        	
	v1[l] = omega*r;
	if(Fluidtype == DUST) {
          v3[l] =  sin(zmin(k))*stokes[id-1]*omega*(r*cos(zmin(k)));
          v2[l] = -cos(zmin(k))*stokes[id-1]*omega*(r*cos(zmin(k)));
        }

	h *= hd;
	if (FLARINGINDEX == 0.0) {
	  rho[l] = sigma/sqrt(2.0*M_PI)/(R0*ASPECTRATIO*hd)*pow(r/R0,-xi)*pow(sin(Zmed(k)),-beta-xi+1./(h*h));
        } else {
          rho[l] = sigma/sqrt(2.0*M_PI)/(R0*ASPECTRATIO*hd)*pow(r/R0,-xi)* \
            pow(sin(Zmed(k)),-xi-beta)*exp((1.-pow(sin(Zmed(k)),-2.*FLARINGINDEX))/2./FLARINGINDEX/(h*h));
        }

	if(Fluidtype == GAS) v1[l] *= sqrt(pow(sin(Zmed(k)),-2.*FLARINGINDEX)-(beta+xi)*h*h);       
	v1[l] -= OMEGAFRAME*r*sin(Zmed(k));


	if(rho[l]<FLOORMIN) rho[l] = FLOORMIN;        	

      }
    }
  }
}



void CondInit() {

  

#ifdef ADIABATIC
      prs_error("WARNING! The gas energy is modified in substep2 and substep3 and therefore it must be communicated to the dust fluids for the DIFFUSION module and DRAGFORCE module when DUSTSIZE is enable. The energy in condinit must be modified accordenly as well.");
#endif
  
  int feedback = YES;
  char dust_name[MAXNAMELENGTH];
  int global_index;
  int id;
  
  
  for (id = 0; id<NFluids_per_rank; id++) {
    global_index = FluidColor*NFluids_per_rank+id;
    if (global_index == 0) {
      Fluids[id] = CreateFluid("gas",GAS);
      SelectFluid(id);
      _CondInit(0);
    }
    else {
      sprintf(dust_name,"dust%d",global_index);
      Fluids[id]  = CreateFluid(dust_name, DUST);
      SelectFluid(id);
      _CondInit(global_index);
    }
  }

}

