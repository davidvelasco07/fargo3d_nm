#include "fargo3d.h"

void FillCPUGrid(tGrid_CPU *desc)
{
  long i, j, k, le, ip[2];
  long idx[3];
  long gncell[3], stride[3];
  real xmin, xmax, ymin, ymax, zmin, zmax;
  real rmin, rmax, phimin, phimax;
  real thetamin, thetamax, dl, coord;
  real volume = 1.0, *InvVolume, **Center, **InterSurface, *Edges[3];
  real *Metric[3][2], *InvMetric[3][2];
  real frac;
  long size, count = 0;
  int global_fluid_id;
  FluidPatch *previousFluid = NULL;
  FluidPatch *fluid;
  Fluid *F3Dfluid;
  real *optical_depth;
  for (i = 0; i < 3; i++)
  {
    gncell[i] = desc->gncell[i];
    stride[i] = desc->stride[i];
  }
#ifdef GPU
  desc->src.nbCommsUp = 0;
  desc->src.nbCommsDown = 0;
  desc->src.nbCommsFlux = 0;
  desc->dst.nbCommsUp = 0;
  desc->dst.nbCommsDown = 0;
  desc->dst.nbCommsFlux = 0;
#endif
  size = gncell[0] * gncell[1] * gncell[2];
  InvVolume = prs_malloc(size * sizeof(real));
  desc->InvVolume = InvVolume;
  InterSurface = multiple_alloc_1D(3, size);
  Center = multiple_alloc_1D(3, size);
  for (i = 0; i < 3; i++)
  { /* One needs to leave 3 here, instead
				   of NDIM, in case one has to refer
				   to a non-indifferent extra
				   dimension if NDIM<3 */
    ip[0] = (i == 0);
    ip[1] = 2 - (i == 2);
    desc->Center[i] = Center[i];
    desc->InterSurface[i] = InterSurface[i];
    Edges[i] = prs_malloc((gncell[i] + 1) * sizeof(real));
    for (le = 0; le < 2; le++)
    {
      Metric[i][le] = prs_malloc(gncell[ip[le]] * sizeof(real));
      desc->Metric[i][le] = Metric[i][le];
      InvMetric[i][le] = prs_malloc(gncell[ip[le]] * sizeof(real));
      desc->InvMetric[i][le] = InvMetric[i][le];
    }
    desc->Edges[i] = Edges[i];
  }
  for (j = 0; j < 3; j++)
  {
    for (i = 0; i <= gncell[j]; i++)
    {
      //frac = (real)(i*desc->Parent->dn[j]+desc->gncorner_min[j])/	\
	//(real)(ncorner_max0[j]-ncorner_min0[j]);
      //Edges[j][i] = frac*(corner_max0[j]-corner_min0[j])+corner_min0[j];
      Edges[j][i] = desc->Parent->Edges[j][i + desc->pcorner_min[j]];
    }
  }
  optical_depth = prs_malloc(gncell[0] * gncell[2] * sizeof(real));
  desc->optical_depth = optical_depth;
  for (k = 0; k < gncell[2]; k++)
  {
    for (j = 0; j < gncell[1]; j++)
    {
      for (i = 0; i < gncell[0]; i++)
      {
        le = i * stride[0] + j * stride[1] + k * stride[2];
        idx[0] = i;
        idx[1] = j;
        idx[2] = k;
#ifdef CARTESIAN
        xmin = Edges[CoordNb[0]][idx[CoordNb[0]]];
        xmax = Edges[CoordNb[0]][idx[CoordNb[0]] + 1];
        ymin = Edges[CoordNb[1]][idx[CoordNb[1]]];
        ymax = Edges[CoordNb[1]][idx[CoordNb[1]] + 1];
        zmin = Edges[CoordNb[2]][idx[CoordNb[2]]];
        zmax = Edges[CoordNb[2]][idx[CoordNb[2]] + 1];
        Center[CoordNb[0]][le] = .5 * (xmin + xmax);
        Center[CoordNb[1]][le] = .5 * (ymin + ymax);
        Center[CoordNb[2]][le] = .5 * (zmin + zmax);
        InterSurface[CoordNb[0]][le] = InterSurface[CoordNb[1]][le] =
            InterSurface[CoordNb[2]][le] = volume = 1.0;
        /* X coordinate */
        dl = (_X_ < NDIM ? xmax - xmin : 1.0);
        volume *= dl;
        InterSurface[CoordNb[1]][le] *= dl;
        InterSurface[CoordNb[2]][le] *= dl;
        /* Y coordinate */
        dl = (_Y_ < NDIM ? ymax - ymin : 1.0);
        volume *= dl;
        InterSurface[CoordNb[0]][le] *= dl;
        InterSurface[CoordNb[2]][le] *= dl;
        /* Z coordinate */
        dl = (_Z_ < NDIM ? zmax - zmin : 1.0);
        volume *= dl;
        InterSurface[CoordNb[0]][le] *= dl;
        InterSurface[CoordNb[1]][le] *= dl;
#endif
#ifdef CYLINDRICAL
        rmin = Edges[CoordNb[0]][idx[CoordNb[0]]];
        rmax = Edges[CoordNb[0]][idx[CoordNb[0]] + 1];
        phimin = Edges[CoordNb[1]][idx[CoordNb[1]]];
        phimax = Edges[CoordNb[1]][idx[CoordNb[1]] + 1];
        zmin = Edges[CoordNb[2]][idx[CoordNb[2]]];
        zmax = Edges[CoordNb[2]][idx[CoordNb[2]] + 1];
        Center[CoordNb[0]][le] = .5 * (rmin + rmax);
        Center[CoordNb[1]][le] = .5 * (phimin + phimax);
        Center[CoordNb[2]][le] = .5 * (zmin + zmax);
        InterSurface[CoordNb[0]][le] = InterSurface[CoordNb[1]][le] =
            InterSurface[CoordNb[2]][le] = volume = 1.0;
        /* Radial coordinate */
        dl = (_RAD_ < NDIM ? .5 * (rmax * rmax - rmin * rmin) : rmin);
        volume *= dl;
        InterSurface[CoordNb[1]][le] *= (_RAD_ < NDIM ? rmax - rmin : 1.0);
        InterSurface[CoordNb[2]][le] *= dl;
        /* Azimuth coordinate */
        dl = (_AZIM_ < NDIM ? phimax - phimin : 1.0);
        volume *= dl;
        InterSurface[CoordNb[0]][le] *= dl * rmin;
        InterSurface[CoordNb[2]][le] *= dl;
        /* Z coordinate */
        dl = (_Z_ < NDIM ? zmax - zmin : 1.0);
        volume *= dl;
        InterSurface[CoordNb[0]][le] *= dl;
        InterSurface[CoordNb[1]][le] *= dl;
#endif
#ifdef SPHERICAL
        rmin = Edges[CoordNb[0]][idx[CoordNb[0]]];
        rmax = Edges[CoordNb[0]][idx[CoordNb[0]] + 1];
        phimin = Edges[CoordNb[1]][idx[CoordNb[1]]];
        phimax = Edges[CoordNb[1]][idx[CoordNb[1]] + 1];
        thetamin = Edges[CoordNb[2]][idx[CoordNb[2]]];
        thetamax = Edges[CoordNb[2]][idx[CoordNb[2]] + 1];
        /* Intersurface is strictly homogoneous to length^2 */
        Center[CoordNb[0]][le] = (2. / 3.) * (rmax * rmax * rmax - rmin * rmin * rmin) / (rmax * rmax - rmin * rmin);
        Center[CoordNb[1]][le] = .5 * (phimin + phimax);
        if (.5 * (thetamin + thetamax) < M_PI / 2)
          Center[CoordNb[2]][le] = asin((cos(thetamin) - cos(thetamax)) / (thetamax - thetamin));
        else
          Center[CoordNb[2]][le] = M_PI - asin((cos(thetamin) - cos(thetamax)) / (thetamax - thetamin));
        InterSurface[CoordNb[0]][le] = InterSurface[CoordNb[1]][le] =
            InterSurface[CoordNb[2]][le] = volume = 1.0;
        /* Radial coordinate */
        dl = (_RAD_ < NDIM ? ONETHIRD * (rmax * rmax * rmax - rmin * rmin * rmin) : rmin * rmin);
        volume *= dl;
        dl = (_RAD_ < NDIM ? .5 * (rmax * rmax - rmin * rmin) : rmin);
        InterSurface[CoordNb[1]][le] *= dl;
        InterSurface[CoordNb[2]][le] *= dl;
        /* Azimuth coordinate */
        dl = (_AZIM_ < NDIM ? phimax - phimin : 1.0);
        volume *= dl;
        InterSurface[CoordNb[0]][le] *= dl * rmin * rmin;
        InterSurface[CoordNb[2]][le] *= dl * sin(thetamin);
        /* colatitude coordinate */
        dl = (_COLAT_ < NDIM ? (cos(thetamin) - cos(thetamax)) : sin(thetamin));
        volume *= dl;
        InterSurface[CoordNb[0]][le] *= dl;
        InterSurface[CoordNb[1]][le] *= (_COLAT_ < NDIM ? thetamax - thetamin : 1.0);
#endif
        InvVolume[le] = 1.0 / volume;
      }
    }
  }
  for (i = 0; i < 3; i++)
  { /* One needs to leave 3 here, instead
				   of NDIM, in case one has to refer
				   to a non-indifferent extra
				   dimension if NDIM<3 */
    ip[0] = (i == 0);
    ip[1] = 2 - (i == 2);
    for (le = 0; le < 2; le++)
    {
      for (j = 0; j < gncell[ip[le]]; j++)
      {
        coord = Center[ip[le]][j * stride[ip[le]]];
#ifdef CARTESIAN
        coord = 1.0;
#endif
#ifdef CYLINDRICAL
        if (!((i == _AZIM_) && (ip[le] == _RAD_)))
          coord = 1.0;
#endif
#ifdef SPHERICAL
        if ((i == _RAD_) || ((i == _COLAT_) && (ip[le] == _AZIM_)))
          coord = 1.0;
        if ((i == _AZIM_) && (ip[le] == _COLAT_))
          coord = sin(coord);
#endif
        Metric[i][le][j] = coord;
        InvMetric[i][le][j] = 1. / coord;
      }
    }
  }

  /* A comment on the Metric[3][2] array that should have been written many years ago... */
  /* In cartesian coordinates this matrix has just ones everywhere.
     In cylindrical coordinates (with COORDPERMUT set to 213) it has the following shape:

    [[r 1]
     [1 1]
     [1 1]]

 and in spherical coordinates (COORDPERMUT set to 213 also):

    [[r sin\theta]
     [1   1      ]
     [1   r      ]] */

  /* This metric matrix gives the length of an arc of a given
coordinate by simply multiplying, on a given line, the first element
by the second. Example: the length of an azimuthal arc of size d\phi
in cylindrical coordinates (azimuth == 1st coord ==> first line of
matrix) is given by: dl = d\phi * r * 1.

Another example: the length of a colatitude arc of size d\theta in
spherical coordinates is given by (3rd line, last coordinate...): 
dl = d\theta * 1 * r.

Another example: the length of an azimuthal arc of size d\phi in
spherical coordinates is given by (1st line, first coordinate...): 
dl = d\phi * r * sin(\theta).

The matrix is correctly evaluated regardless of the value of COORDPERMUT.
It should always feature in pair, involving the multiplication of the first
element by the second element. So we should always see a [0] ... [1] ...
in lines involving the metric.

The reason why we have to multiply the elements (instead of storing a single
array that is the result of the multiplication) is that we only need, this
way, to store 1D array, with a much lighter memory imprint.

The array InvMetric has just coefficients which are 1 by 1 the inverse of those
of Metric.
  */
  tGrid_CPU *grid = desc;
  grid->Xmed = prs_malloc(sizeof(real) * grid->gncell[_X_]);
  grid->Ymed = prs_malloc(sizeof(real) * grid->gncell[_Y_]);
  grid->Zmed = prs_malloc(sizeof(real) * grid->gncell[_Z_]);

  grid->InvDiffXmed = prs_malloc(sizeof(real) * (grid->gncell[_X_] + 1));
  grid->InvDiffYmed = prs_malloc(sizeof(real) * (grid->gncell[_Y_] + 1));
  grid->InvDiffZmed = prs_malloc(sizeof(real) * (grid->gncell[_Z_] + 1));

  grid->Sxj = prs_malloc(sizeof(real) * grid->gncell[_Y_]);
  grid->Syj = prs_malloc(sizeof(real) * grid->gncell[_Y_]);
  grid->Szj = prs_malloc(sizeof(real) * grid->gncell[_Y_]);

  grid->Sxk = prs_malloc(sizeof(real) * grid->gncell[_Z_]);
  grid->Syk = prs_malloc(sizeof(real) * grid->gncell[_Z_]);
  grid->Szk = prs_malloc(sizeof(real) * grid->gncell[_Z_]);

  grid->InvVj = prs_malloc(sizeof(real) * grid->gncell[_Y_]);

  Xmin = grid->Edges[_X_];
  Ymin = grid->Edges[_Y_];
  Zmin = grid->Edges[_Z_];
  Xmed = grid->Xmed;
  Ymed = grid->Ymed;
  Zmed = grid->Zmed;
  InvDiffXmed = grid->InvDiffXmed;
  InvDiffYmed = grid->InvDiffYmed;
  InvDiffZmed = grid->InvDiffZmed;
  Sxj = grid->Sxj;
  Syj = grid->Syj;
  Szj = grid->Szj;
  Sxk = grid->Sxk;
  Syk = grid->Syk;
  Szk = grid->Szk;
  InvVj = grid->InvVj;

  Nx = grid->ncell[_X_];
  Ny = grid->ncell[_Y_];
  Nz = grid->ncell[_Z_];

  Dx = (Xmin[NGHX + Nx] - Xmin[NGHX]) / Nx;

  for (i = 0; i < Nx + 2 * NGHX; i++)
  {
    Xmed(i) = 0.5 * (Xmin(i + 1) + Xmin(i));
  }
  for (j = 0; j < Ny + 2 * NGHY; j++)
  {
    Ymed(j) = 0.5 * (Ymin(j + 1) + Ymin(j));
  }
  for (k = 0; k<Nz+2*NGHZ; k++) {
#ifdef Z
    Zmed(k) = 0.5*(Zmin(k+1)+Zmin(k));
#else
    Zmed(k) = Zmin(k+1) = Zmin(k) = 0;
#endif
  }

  InitSurfaces();

#ifdef GPU
  DevMalloc(&grid->Xmin_d, sizeof(real) * (Nx + 2 * NGHX + 1));
  DevMemcpyH2D(grid->Xmin_d, Xmin, sizeof(real) * (Nx + 2 * NGHX + 1));
  DevMalloc(&grid->Ymin_d, sizeof(real) * (Ny + 2 * NGHY + 1));
  DevMemcpyH2D(grid->Ymin_d, Ymin, sizeof(real) * (Ny + 2 * NGHY + 1));
  DevMalloc(&grid->Zmin_d, sizeof(real) * (Nz + 2 * NGHZ + 1));
  DevMemcpyH2D(grid->Zmin_d, Zmin, sizeof(real) * (Nz + 2 * NGHZ + 1));

  DevMalloc(&grid->Sxj_d, sizeof(real) * (Ny + 2 * NGHY));
  DevMemcpyH2D(grid->Sxj_d, Sxj, sizeof(real) * (Ny + 2 * NGHY));
  DevMalloc(&grid->Syj_d, sizeof(real) * (Ny + 2 * NGHY));
  DevMemcpyH2D(grid->Syj_d, Syj, sizeof(real) * (Ny + 2 * NGHY));
  DevMalloc(&grid->Szj_d, sizeof(real) * (Ny + 2 * NGHY));
  DevMemcpyH2D(grid->Szj_d, Szj, sizeof(real) * (Ny + 2 * NGHY));

  DevMalloc(&grid->Sxk_d, sizeof(real) * (Nz + 2 * NGHZ));
  DevMemcpyH2D(grid->Sxk_d, Sxk, sizeof(real) * (Nz + 2 * NGHZ));
  DevMalloc(&grid->Syk_d, sizeof(real) * (Nz + 2 * NGHZ));
  DevMemcpyH2D(grid->Syk_d, Syk, sizeof(real) * (Nz + 2 * NGHZ));
  DevMalloc(&grid->Szk_d, sizeof(real) * (Nz + 2 * NGHZ));
  DevMemcpyH2D(grid->Szk_d, Szk, sizeof(real) * (Nz + 2 * NGHZ));

  DevMalloc(&grid->InvVj_d, sizeof(real) * (Ny + 2 * NGHY));
  DevMemcpyH2D(grid->InvVj_d, InvVj, sizeof(real) * (Ny + 2 * NGHY));
#endif
  char dust_name[MAXNAMELENGTH];
  for (i = 0; i < NFluids_per_rank; i++)
  {
    global_fluid_id = FluidColor * NFluids_per_rank + i;
    if (global_fluid_id == 0)
      desc->fluid = CreateFluidPatch(desc, "gas", GAS);
    else
    {
      sprintf(dust_name, "dust%d", global_fluid_id);
      desc->fluid = CreateFluidPatch(desc, dust_name, DUST);
    }
    SelectGrid(desc);
    fluid = desc->fluid;
    F3Dfluid = CreateFluid(fluid->Name, fluid->Fluidtype);
    F3Dfluid->Density->field_cpu = fluid->Density->Field;
    F3Dfluid->Energy->field_cpu = fluid->Energy->Field;
#ifdef X
    F3Dfluid->Vx->field_cpu = fluid->Velocity->Field[_X_];
    F3Dfluid->Vx_temp->field_cpu = fluid->V_temp->Field[_X_];
#endif
#ifdef Y
    F3Dfluid->Vy->field_cpu = fluid->Velocity->Field[_Y_];
    F3Dfluid->Vy_temp->field_cpu = fluid->V_temp->Field[_Y_];
#endif
#ifdef Z
    F3Dfluid->Vz->field_cpu = fluid->Velocity->Field[_Z_];
    F3Dfluid->Vz_temp->field_cpu = fluid->V_temp->Field[_Z_];
#endif
    for (j = 0; j < 3; j++)
      for (k = 0; k < 2; k++)
        F3Dfluid->Fluxes[j][k] = fluid->Fluxes[j][k];
#ifdef STOCKHOLM
    F3Dfluid->Density0->field_cpu = fluid->Rho0->Field;
#ifdef X
    F3Dfluid->Vx0->field_cpu = fluid->Vx0->Field;
#endif
#ifdef Y
    F3Dfluid->Vy0->field_cpu = fluid->Vy0->Field;
#endif
#ifdef Z
    F3Dfluid->Vz0->field_cpu = fluid->Vz0->Field;
#endif
#endif

#ifdef GPU
    F3Dfluid->Density->field_gpu = fluid->Density->Field_gpu;
    F3Dfluid->Density->fresh_cpu = &(fluid->Density->fresh_cpu);
    F3Dfluid->Density->fresh_gpu = &(fluid->Density->fresh_gpu);
    
#ifdef STOCKHOLM
    F3Dfluid->Density0->field_gpu = fluid->Rho0->Field_gpu;
    F3Dfluid->Density0->fresh_cpu = &(fluid->Rho0->fresh_cpu);
    F3Dfluid->Density0->fresh_gpu = &(fluid->Rho0->fresh_gpu);
#endif
#ifdef X
    F3Dfluid->Vx->field_gpu = fluid->Velocity->Field_gpu[_X_];
    F3Dfluid->Vx->fresh_cpu = &(fluid->Velocity->fresh_cpu[_X_]);
    F3Dfluid->Vx->fresh_gpu = &(fluid->Velocity->fresh_gpu[_X_]);
    F3Dfluid->Vx_temp->field_gpu = fluid->V_temp->Field_gpu[_X_];
    F3Dfluid->Vx_temp->fresh_cpu = &(fluid->V_temp->fresh_cpu[_X_]);
    F3Dfluid->Vx_temp->fresh_gpu = &(fluid->V_temp->fresh_gpu[_X_]);
#ifdef STOCKHOLM
    F3Dfluid->Vx0->field_gpu = fluid->Vx0->Field_gpu;
    F3Dfluid->Vx0->fresh_cpu = &(fluid->Vx0->fresh_cpu);
    F3Dfluid->Vx0->fresh_gpu = &(fluid->Vx0->fresh_gpu);
#endif
#endif

#ifdef Y
    F3Dfluid->Vy->field_gpu = fluid->Velocity->Field_gpu[_Y_];
    F3Dfluid->Vy->fresh_cpu = &(fluid->Velocity->fresh_cpu[_Y_]);
    F3Dfluid->Vy->fresh_gpu = &(fluid->Velocity->fresh_gpu[_Y_]);
    F3Dfluid->Vy_temp->field_gpu = fluid->V_temp->Field_gpu[_Y_];
    F3Dfluid->Vy_temp->fresh_cpu = &(fluid->V_temp->fresh_cpu[_Y_]);
    F3Dfluid->Vy_temp->fresh_gpu = &(fluid->V_temp->fresh_gpu[_Y_]);
#ifdef STOCKHOLM
    F3Dfluid->Vy0->field_gpu = fluid->Vy0->Field_gpu;
    F3Dfluid->Vy0->fresh_cpu = &(fluid->Vy0->fresh_cpu);
    F3Dfluid->Vy0->fresh_gpu = &(fluid->Vy0->fresh_gpu);
#endif
#endif

#ifdef Z
    F3Dfluid->Vz->field_gpu = fluid->Velocity->Field_gpu[_Z_];
    F3Dfluid->Vz->fresh_cpu = &(fluid->Velocity->fresh_cpu[_Z_]);
    F3Dfluid->Vz->fresh_gpu = &(fluid->Velocity->fresh_gpu[_Z_]);
    F3Dfluid->Vz_temp->field_gpu = fluid->V_temp->Field_gpu[_Z_];
    F3Dfluid->Vz_temp->fresh_cpu = &(fluid->V_temp->fresh_cpu[_Z_]);
    F3Dfluid->Vz_temp->fresh_gpu = &(fluid->V_temp->fresh_gpu[_Z_]);
#ifdef STOCKHOLM
    F3Dfluid->Vz0->field_gpu = fluid->Vz0->Field_gpu;
    F3Dfluid->Vz0->fresh_cpu = &(fluid->Vz0->fresh_cpu);
    F3Dfluid->Vz0->fresh_gpu = &(fluid->Vz0->fresh_gpu);
#endif
#endif
    F3Dfluid->Energy->field_gpu = fluid->Energy->Field_gpu;
    F3Dfluid->Energy->fresh_cpu = &(fluid->Energy->fresh_cpu);
    F3Dfluid->Energy->fresh_gpu = &(fluid->Energy->fresh_gpu);
#ifdef STOCKHOLM
    F3Dfluid->Energy0->field_gpu = fluid->Energy0->Field_gpu;
    F3Dfluid->Energy0->fresh_cpu = &(fluid->Energy0->fresh_cpu);
    F3Dfluid->Energy0->fresh_gpu = &(fluid->Energy0->fresh_gpu);
#endif

#endif
    desc->Fluids[i] = F3Dfluid;
    fluid->FluidRank = global_fluid_id;
    F3Dfluid->FluidRank = global_fluid_id;
    desc->fluid->next = previousFluid;
    if (previousFluid != NULL)
      previousFluid->prev = desc->fluid;
    previousFluid = desc->fluid;
  }
}
