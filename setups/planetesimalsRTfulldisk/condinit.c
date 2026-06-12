#include "fargo3d.h"

void CondInit() {
  int i,j,k;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *er;
  real *rho;
  real h;
  
  real omega;
  real r, r3, temperature;
  real rplanet;
  
  KAPPA	       *= (R0*R0/MSTAR)/(R0_CGS*R0_CGS/MSTAR_CGS);
  DOUBLINGTIME *= ((sqrt(R0*R0*R0/G/MSTAR))/(sqrt(R0_CGS*R0_CGS*R0_CGS/G_CGS/MSTAR_CGS)))*(1.0e6*365*86400);
  RHOSOLID *= (MSTAR/(R0*R0*R0))/(MSTAR_CGS/(R0_CGS*R0_CGS*R0_CGS));

  rplanet = pow(PLANETMASS*3.0/(4.0*M_PI)/RHOSOLID,0.33333333333);
  HeatingRate = G*PLANETMASS*PLANETMASS/DOUBLINGTIME/rplanet;

  if (CPU_Master) {
    printf("\nACCRETION PROPERTIES:\n");
    printf("=====================\n");
    printf("PLANET RADIUS:\t%e Earth radii\n", rplanet/R0*R0_CGS/6371e5);
    printf("LUMINOSITY:\t%e erg/s\n", HeatingRate/(pow(G*MSTAR/R0,2.5)/G)*(pow(G_CGS*MSTAR_CGS/R0_CGS,2.5)/G_CGS));
    printf("           \t%e Solar Units\n", HeatingRate/(pow(G*MSTAR/R0,2.5)/G)*(pow(G_CGS*MSTAR_CGS/R0_CGS,2.5)/G_CGS)/3.826e33);
    printf("           \t%3.1e Jupiter Units\n", HeatingRate/(pow(G*MSTAR/R0,2.5)/G)*(pow(G_CGS*MSTAR_CGS/R0_CGS,2.5)/G_CGS)/3.826e33/8.7e-10);


    printf("\n");
  }

  INSPECT_REAL(DT);
  INSPECT_REAL(PLANETMASS);
  INSPECT_REAL(rplanet);
  INSPECT_REAL(RHOSOLID);
  INSPECT_REAL(KAPPA);
  INSPECT_REAL(HeatingRate);
  INSPECT_REAL(DOUBLINGTIME);

  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  er  = Energyrad->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;

  real aspectratio;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {

      aspectratio = pow(9.0/32.0*SIGMA0*SIGMA0*KAPPA*NU*pow(R_MU,4)/STEFANK/pow(G*MSTAR,3)*R0,0.125)*pow(Ymed(j)/R0,(1.0-2*SIGMASLOPE)/8.0)*1.1;
      

      h = aspectratio*Ymed(j);

      r = Ymed(j);
      r3 = r*r*r;
      omega = sqrt(G*MSTAR/(r3));
      for (i=0; i<Nx; i++) {
	v2[l] = v3[l] = 0.0;
	v1[l] = omega*r;
#ifdef CYLINDRICAL
	rho[l] = SIGMA0*pow(r/R0,-SIGMASLOPE)*exp(-pow(Zmed(k)/h,2.0)/2.0)/(ZMAX-ZMIN);
#else
	real xi = SIGMASLOPE+1.+FLARINGINDEX;
	real beta = 1.-2*FLARINGINDEX;
	real h = aspectratio*pow(r/R0,FLARINGINDEX);
	if (FLARINGINDEX == 0.0) {
	  rho[l] = SIGMA0/sqrt(2.0*M_PI)/(R0*aspectratio)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-beta-xi+1./(h*h));
	} else {
	  rho[l] = SIGMA0/sqrt(2.0*M_PI)/(R0*aspectratio)*pow(r/R0,-xi)* \
	    pow(sin(Zmed(k)),-xi-beta)*					\
	    exp((1.-pow(sin(Zmed(k)),-2.*FLARINGINDEX))/2./FLARINGINDEX/(h*h));
	}
	v1[l] *= sqrt(pow(sin(Zmed(k)),-2.*FLARINGINDEX)-(beta+xi)*h*h);
	v1[l] -= OMEGAFRAME*r*sin(Zmed(k));
#endif
	e[l] = rho[l]*h*h*G*MSTAR/r/(GAMMA-1.0);
	temperature = (GAMMA-1.0)*e[l]/(rho[l]*R_MU);
	er[l] = 4.0*STEFANK*pow(temperature,4.0)/C0; // Initial estimate of radiative energy
      }
    }
  }
}
