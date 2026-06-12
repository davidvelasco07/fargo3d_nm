#include "fargo3d.h"

void InitDensPlanet() {

  int i,j,k;
  real *field;

  field = Density->field_cpu;
    
  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = 0;
  int end_i = Nx;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
	field[l] = SIGMA0*pow((Ymed(j)/R0),-SIGMASLOPE)*(1.0+NOISE*(drand48()-.5));
      }
    }
  }
}

void InitSoundSpeedPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t, omega, vk;
  real rho;
  FILE *fo;
  real *d;
  real *e;

  field = Energy->field_cpu;
  d = Density->field_cpu;

  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = 0;
  int end_i = Nx;
  
  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {	
	r = Ymed(j);
	vk = sqrt(G*MSTAR/r);
	field[l] = ASPECTRATIO * pow(Ymed(j)/R0, FLARINGINDEX) * vk; //sqrt(G*MSTAR/Ymed(j))
#ifdef ADIABATIC
	field[l] = field[l]*field[l]*d[l]/(GAMMA-1.0);
#endif
      }
    }
  }    
}

void CentrifugalBalance() {

  INPUT(Pot);
  INPUT(Density);
  INPUT(Pressure);

  printf("WARNING: Forcing to machine precision the radial forces balance.\n");
  printf("================================================================\n\n");

  int i,j,k;
  real grad;

  real* pot = Pot->field_cpu;
  real* p = Pressure->field_cpu;
  real* rho = Density->field_cpu;
  real* vphi = Vx->field_cpu;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=1; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx; i++) {
	grad = 2.0*(p[l]-p[lym])/(rho[l]+rho[lym])/(ymed(j)-ymed(j-1)) + (pot[l]-pot[lym])/(ymed(j)-ymed(j-1));
	vphi[l] = 2.0*(sqrt(ymin(j)*grad) - ymin(j)*OMEGAFRAME) - vphi[lym];
      }
    }
  }
}

void InitVazimPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, t;
  FILE *fo;
  real vt, omega;
  real *vr;
  real *cs;
  real *p;
  real *pot;
  real *rho;

  field = Vx->field_cpu;
  vr = Vy->field_cpu;
  cs = Energy->field_cpu;
  pot = Pot->field_cpu;
  p = Pressure->field_cpu;
  rho = Density->field_cpu;

  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = 0;
  int end_i = Nx;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {
	r = Ymed(j);
	omega = sqrt(G*MSTAR/r/r/r);
	vt = omega*r*sqrt(1.0+pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
			  (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));	
	vt -= OMEGAFRAME*r;
	field[l] = vt*(1.+ASPECTRATIO*NOISE*(drand48()-.5));
	vr[l] = r*omega*ASPECTRATIO*NOISE*(drand48()-.5);	
      }
    }
  }
}

void CondInit() {
  
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);

  int i,j,k;
  int index;
  real vt;
  real *field;
  real *rho;
  real *v1;
  real *v2;
  real *e;

#ifdef PLANETS
  Sys = InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();
  if(COROTATING)
    OMEGAFRAME = GetPsysInfo(FREQUENCY);
  else
#endif
    OMEGAFRAME = OMEGAFRAME;
  
  InitDensPlanet ();
  InitSoundSpeedPlanet ();  
  InitVazimPlanet ();
  
  if (CENTRIFUGALBALANCE){
    Potential_cpu();
#ifdef ADIABATIC
    FARGO_SAFE(ComputePressureFieldAd_cpu());
#endif
    
#ifdef ISOTHERMAL
    FARGO_SAFE(ComputePressureFieldIso_cpu());
#endif
    CentrifugalBalance();
  }

}
