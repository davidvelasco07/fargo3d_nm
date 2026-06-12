#include "fargo3d.h"

void CentrifugalBalance() {

  INPUT(Pot);
  INPUT(Density);
  INPUT(Pressure);

  printf("WARNING: Forcing to machine precision the radial balance.\n");
  printf("================================================================\n\n");

  int i,j,k;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;

  real grad;

  real* pot = Pot->field_cpu;
  real* p = Pressure->field_cpu;
  real* rho = Density->field_cpu;
  real* vphi = Vx->field_cpu;

#ifdef MHD
  real* bx = Bx->field_cpu;
  real* bz = Bz->field_cpu;
  real db1;
  real db2;
  real bmeanm;
  real bmean;
#endif

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=1; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx; i++) {
	grad = 2.0*(p[l]-p[lym])/(rho[l]+rho[lym])/(ymed(j)-ymed(j-1)) + (pot[l]-pot[lym])/(ymed(j)-ymed(j-1));
#ifdef MHD
	//MAGNETIG PRESSURE
	bmean  = 0.5*(bx[l] + bx[lxp]);
	bmeanm = 0.5*(bx[lym] + bx[lxp-pitch]);
	db1 = (bmean*bmean-bmeanm*bmeanm);
	bmean  = 0.5*(bz[l] + bz[lzp]);
	bmeanm = 0.5*(bz[lym] + bz[lym+stride]);
	db2 = (bmean*bmean-bmeanm*bmeanm);
	grad += 1.0*(db1 + db2)/(rho[l]+rho[lym])/(ymed(j)-ymed(j-1))/MU0;
#endif
	vphi[l] = 2.0*(sqrt(ymin(j)*grad) - ymin(j)*OMEGAFRAME) - vphi[lym];
      }
    }
  }
}

void CondInit() {

  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);
#ifdef MHD
  OUTPUT(Bx);
  OUTPUT(By);
#endif

  srand48 (NOISESEED); //Same seed for all processes

  int i,j,k;
  real vphi_t;
  real r,omega,h;
  real b0;

  real *rho   = Density->field_cpu;
  real *cs    = Energy->field_cpu;
  real *vphi  = Vx->field_cpu;
  real *vr    = Vy->field_cpu;
  real *vz    = Vz->field_cpu;
  real *bx    = Bx->field_cpu;
  real *by    = By->field_cpu;
  real *bz    = Bz->field_cpu;

  real mdot;

  mdot = (SIGMA0/(ZMAX-ZMIN))*3*M_PI*5e-2*ASPECTRATIO*ASPECTRATIO;
  printf("%--------------------------------------------------------------lf\n",mdot);
  rho = Density->field_cpu;
  //  b0 = sqrt(2.0*SIGMA0*ASPECTRATIO*MU0/(BETA*sqrt(2.0*M_PI)));
  b0 = ASPECTRATIO*sqrt(2.0*SIGMA0*MU0/BETA);

  for (k = 0; k<Nz+2*NGHZ; k++) {
    for (j = 0; j<Ny+2*NGHY; j++) {
      for (i = 0; i<Nx; i++) {
	r = Ymed(j);
	omega = sqrt(G*MSTAR/r/r/r);
	rho[l] = SIGMA0/(ZMAX-ZMIN)*pow((r/R0),-SIGMASLOPE);
	cs[l] = ASPECTRATIO * pow(r/R0, FLARINGINDEX) * sqrt(G*MSTAR/r);
#ifdef ADIABATIC
	cs[l] = rho[l]*rho[l]*d[l]/(GAMMA-1.0);
#endif
	h = cs[l]/(omega*r);

	vphi_t  = omega*r*sqrt(1.0+(-SIGMASLOPE + (2.0*FLARINGINDEX-1.0) + 2.0*(MAGNETICSLOPE+1)/BETA)*h*h);
	vphi_t -= OMEGAFRAME*r;
	vphi[l] = vphi_t*(1.+ASPECTRATIO*NOISE*(drand48()-.5));
	vr[l]   = r*omega*ASPECTRATIO*NOISE*(drand48()-.5);
	
	vr[l]   += -mdot/(2*M_PI*Ymin(j)*(SIGMA0/(ZMAX-ZMIN)*pow((Ymin(j)/R0),-SIGMASLOPE)));
	//	printf("%lf %lf\n",vphi[l],vr[l]);
	vz[l]   = r*omega*ASPECTRATIO*NOISE*(drand48()-.5);
	//	bx[l] = cs[l]*sqrt(2.0*rho[l]*MU0/BETA);
	bx[l] = b0*pow((r/R0),MAGNETICSLOPE)*(1.0+NOISE*(drand48()-.5));
	by[l] = 0.0;
	bz[l] = 0.0;
      }
    }
  }
  //  exit(1);
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
