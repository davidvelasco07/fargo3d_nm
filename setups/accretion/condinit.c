#include "fargo3d.h"

void _CondInit(int id)
{
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);
  OUTPUT(Vz);

  int i, j, k, n;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *rho;

  real omega, rho_s, alp;
  real r, r3, xi, beta, h;
  real latitud_m, latitud_p, h_local, efm, efp;

  rho = Density->field_cpu;
  e = Energy->field_cpu;
  v1 = Vx->field_cpu;
  v2 = Vy->field_cpu;
  v3 = Vz->field_cpu;

  real stokes_plus[NFLUIDS];
  real stokes[NFLUIDS - 1];
  real epsilons[NFLUIDS - 1];
  real sq = SQ;
  real slope, hd, sigma;
  // Stokes numbers
  real smax = SMAX;
  real smin = SMIN;
  real ds = (log(smax) - log(smin)) / (NFLUIDS - 1);
  for (n = 0; n < NFLUIDS; n++)
  {
    stokes_plus[n] = smin * exp(ds * n);
  }

  // Density distribution
  slope = 4 - sq;
  for (n = 0; n < NFLUIDS - 1; n++)
  {
    if (slope != 0.)
    {
      epsilons[n] = pow(stokes_plus[n + 1], slope) - pow(stokes_plus[n], slope);
      epsilons[n] *= EPSILON / (pow(smax, slope) - pow(smin, slope));
    }
    else
    {
      epsilons[n] = log(stokes_plus[n + 1] / stokes_plus[n]);
      epsilons[n] *= EPSILON / log(smax / smin);
    }
    stokes[n] = stokes_plus[n + 1];
    if (NFLUIDS == 2)
      stokes[n] = SMAX;
  }
  
  sigma = SIGMA0;
  if (id > 0)
  {
    sigma = SIGMA0 * epsilons[id - 1];
#ifdef STOKESNUMBER
    Coeffval[0] = 1.0 / stokes[id - 1];
#endif

  }

  if (Fluidtype == DUST)
    hd = MIN(1.0, sqrt((DELTA + NU) / (stokes[id - 1] + (DELTA + NU))));
  if (Fluidtype == GAS)
    hd = 1.0;

  xi = SIGMASLOPE + 1. + FLARINGINDEX;
  beta = 1. - 2 * FLARINGINDEX;

  for (k = 0; k < Nz + 2 * NGHZ; k++)
  {
    for (j = 0; j < Ny + 2 * NGHY; j++)
    {
      r = Ymed(j);
      r3 = r * r * r;
      omega = sqrt(G * MSTAR / (r3));

      h = ASPECTRATIO * pow(r / R0, FLARINGINDEX);

      for (i = 0; i < Nx + 2 * NGHX; i++)
      {
        v2[l] = v3[l] = 0.0;
        v1[l] = omega * r;
        if (Fluidtype == DUST)
        {
          v3[l] = sin(zmin(k)) * stokes[id - 1] * omega * (r * cos(zmin(k))) * 0;
          v2[l] = -cos(zmin(k)) * stokes[id - 1] * omega * (r * cos(zmin(k))) * 0;
        }

        if (FLARINGINDEX == 0.0)
        {
          rho[l] = sigma / sqrt(2.0 * M_PI) / (R0 * ASPECTRATIO * hd) * pow(r / R0, -xi) * pow(sin(Zmed(k)), -beta - xi + 1. / (h * hd * h * hd));
        }
        else
        {
          rho[l] = sigma / sqrt(2.0 * M_PI) / (R0 * ASPECTRATIO * hd) * pow(r / R0, -xi) *
                   pow(sin(Zmed(k)), -xi - beta) * exp((1. - pow(sin(Zmed(k)), -2. * FLARINGINDEX)) / 2. / FLARINGINDEX / (h * hd * h * hd));
        }
	if (Fluidtype == DUST) { /* We overwrite the density with an implementation that will correctly handle vertically unresolved discs */
	  h_local = ASPECTRATIO * pow(Ymed(j)/R0, FLARINGINDEX);
	  h_local *= hd;
	  latitud_m = M_PI/2.-Zmin[k+1];
	  latitud_p = M_PI/2.-Zmin[k];
	  efm = erf(latitud_m/h_local/sqrt(2.));
	  efp = erf(latitud_p/h_local/sqrt(2.));
	  rho[l] = sigma*pow(Ymed(j)/R0,-SIGMASLOPE) * (efp-efm) / zone_size_z(j,k) / 2.0;
	}

        if (Fluidtype == GAS)
          v1[l] *= sqrt(pow(sin(Zmed(k)), -2. * FLARINGINDEX) - (beta + xi) * h * h);

        v1[l] -= OMEGAFRAME * r * sin(Zmed(k));

        if (rho[l] < FLOORMIN)
          rho[l] = FLOORMIN;
        
        if(Fluidtype == GAS){
        #ifdef ISOTHERMAL
	        e[l] = h*sqrt(G*MSTAR/r);
        #else
	        e[l] = rho[l]*h*h*G*MSTAR/r/(GAMMA-1.0);
        #endif
	      }
	      else{
	        e[l] = 0.;
	      }
      }
    }
  }
}

void CondInit()
{
  int feedback = YES;
  int id;
  for (id = 0; id < NFluids_per_rank; id++)
  {
    SelectFluid(id);
    _CondInit(Current_Fluid);
  }
}
