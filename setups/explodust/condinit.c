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

void InitVazimPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t;
  real rho;
  FILE *fo;
  real vt, omega;
  real *vr;
  real *cs;

  field = Vx->field_cpu;
  vr = Vy->field_cpu;
  cs = Energy->field_cpu;
    
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
  real qoverhcube;
  FILE *diag;
  char diagfile[2048];

  int index_mass;
  int index_h;
  int SquareDim = 15;

  index_mass = ArrayNb / SquareDim;
  index_h    = ArrayNb % SquareDim;

  

#ifdef PLANETS
  Sys = InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();
  if(COROTATING)
    OMEGAFRAME = GetPsysInfo(FREQUENCY);
  else
#endif
    OMEGAFRAME = OMEGAFRAME;

  ASPECTRATIO = (0.03+(0.10-0.04)*(double)index_h/((double)SquareDim-1.0));
  qoverhcube = exp(log(0.1)+(log(5./0.1))*(double)index_mass/((double)SquareDim-1.0));
  Sys->mass[0] = qoverhcube*ASPECTRATIO*ASPECTRATIO*ASPECTRATIO;
  NU = ALPHAVISCOSITY*ASPECTRATIO*ASPECTRATIO;
 
  sprintf (diagfile, "%s/localparams.log", OUTPUTDIR);
  diag = fopen (diagfile, "w");
  if (diag == NULL) {
    fprintf (stderr, "Could not open localparams.log file\n");
    prs_exit (1);
  }

  fprintf (diag, "Run #%d\n", ArrayNb);
  fprintf (diag, "Index mass: %d\n", index_mass);
  fprintf (diag, "Index h: %d\n", index_h);
  fprintf (diag, "Aspect ratio: %g\n", ASPECTRATIO);
  fprintf (diag, "q over h cube: %g\n", qoverhcube);
  fprintf (diag, "Planet mass: %g\n",  Sys->mass[0]);
  fprintf (diag, "Viscosity: %g\n", NU);
  fclose (diag);

  InitDensPlanet ();
  InitSoundSpeedPlanet ();
  InitVazimPlanet ();
}
