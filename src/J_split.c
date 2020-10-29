#include "fargo3d.h"

static boolean isprime[MAXPRIME];
static long prime[MAXPRIME];
long cpugrid_number;
static boolean primedone = FALSE;

boolean Buildprime()
{
  long i, j;
  for (i = 0; i < MAXPRIME; i++)
  {
    isprime[i] = TRUE;
  }
  for (i = 2; i < MAXPRIME; i++)
  {
    if (isprime[i])
    {
      for (j = 2 * i; j < MAXPRIME; j++)
      {
        if (!(j % i))
          isprime[j] = FALSE;
      }
    }
  }
  isprime[0] = isprime[1] = j = 0;
  for (i = 2; i < MAXPRIME; i++)
  {
    if (isprime[i])
      prime[j++] = i;
  }
  return TRUE;
}

void Primefactors(long n, long *factors, long *nfact)
{
  long j = 0, k = 0;
  if (!primedone)
    primedone = Buildprime();
  while (n > 1)
  {
    if (n % prime[j])
    {
      j++;
    }
    else
    {
      n /= prime[j];
      factors[k++] = prime[j];
    }
  }
  *nfact = k;
}

void Repartition(long *nx, long ncpu, long *MX, long level)
{
  long factors[MAXPRIME], nfact, i, ncomb, pow3[MAXPRIME], j, comb[MAXPRIME], nb;
  long mx[3], p[3], idx;
  real best = 1E30;
  Primefactors(ncpu, factors, &nfact);
  for (i = 0; i <= nfact; i++)
  {
    pow3[i] = (long)(pow(3.0, (real)i) + .01);
  }
  ncomb = pow3[nfact];
  for (i = 0; i < ncomb; i++)
  {
    nb = i;
    for (j = nfact - 1; j >= 0; j--)
    {
      comb[j] = nb / pow3[j];
      nb -= comb[j] * pow3[j];
    }
    for (j = 0; j < 3; j++)
      mx[j] = 1;
    for (j = 0; j < nfact; j++)
      mx[comb[j]] *= factors[j];
    for (j = 0; j < 3; j++)
    {
      p[j] = nx[j] / mx[j];
    }
    idx = labs(p[0] - p[1]) + labs(p[0] - p[2]) + labs(p[1] - p[2]);
#ifndef STANDARD
    if ((level == 0) && (mx[0] > 1))
      idx = best * 2;
#endif
    if (idx < best)
    {
      best = idx;
      for (j = 0; j < 3; j++)
        MX[j] = mx[j];
    }
  }
  ////Hardwiring parallelization on Y direction
  ////This has to be fixed
  for (j = 0; j < 3; j++)
        MX[j] = 1;
  MX[1] = ncpu;
}

void splitgrid(tGrid *grid)
{
  long ncpus, mx[3], i[3], cpu, ii, j, resol_ratio[3];
  long euclid_ratio[3], euclid_rest[3];
  tGrid_CPU *gridcpu;
  ncpus = NDomains;
  Repartition(grid->ncell, ncpus, mx, grid->level);
  for (j = 0; j < 3; j++)
  {
    /* The resol_ratio variable has been introduced so that the interface
       between two CPUgrids at level l be also an interface of level l-1
       (otherwise the ghost buffer is mangled) */
    resol_ratio[j] = (Refine[j] && (grid->level > 0) ? 2L : 1L);
    euclid_ratio[j] = (grid->ncell[j] / resol_ratio[j]) / mx[j];
    euclid_rest[j] = (grid->ncell[j] / resol_ratio[j]) - euclid_ratio[j] * mx[j];
    grid->Ncpus[j] = mx[j];
  }
  //printf("Grid %d, CPUS (%d,%d,%d)\n",grid->number,mx[0],mx[1],mx[2]);
  for (i[0] = 0; i[0] < mx[0]; i[0]++)
  {
    for (i[1] = 0; i[1] < mx[1]; i[1]++)
    {
      for (i[2] = 0; i[2] < mx[2]; i[2]++)
      {
        gridcpu = prs_malloc(sizeof(tGrid_CPU));
        cpu = i[0] + i[1] * mx[0] + i[2] * mx[0] * mx[1];
        gridcpu->color = i[0] + i[2] * mx[0]; //Unique radial beam of CPUs identifier
        gridcpu->key = i[1];
        gridcpu->colorz = i[2];
        gridcpu->parent = grid->number;
        gridcpu->Parent = grid;
        gridcpu->number = cpugrid_number++;
        gridcpu->level = grid->level;
        gridcpu->cpu = cpu;
        for (j = 0; j < 3; j++)
        {
          gridcpu->ncell[j] = euclid_ratio[j];
          if (i[j] < euclid_rest[j])
            gridcpu->ncell[j]++;
          gridcpu->ncell[j] *= resol_ratio[j];
          gridcpu->gncell[j] = gridcpu->ncell[j] + 2 * Nghost[j];
          gridcpu->nsize[j] = gridcpu->ncell[j] * grid->dn[j];
          gridcpu->pcorner_min[j] = i[j] * euclid_ratio[j];
          gridcpu->pcorner_min[j] += (i[j] < euclid_rest[j] ? i[j] : euclid_rest[j]);
          gridcpu->pcorner_min[j] *= resol_ratio[j];
          gridcpu->pcorner_max[j] = gridcpu->pcorner_min[j] + gridcpu->ncell[j];
          gridcpu->corner_min[j] = grid->Edges[j][Nghost[j] + gridcpu->pcorner_min[j]];
          gridcpu->corner_max[j] = grid->Edges[j][Nghost[j] + gridcpu->pcorner_max[j]];
          gridcpu->ncorner_min[j] = gridcpu->pcorner_min[j] * grid->dn[j] + grid->ncorner_min[j];
          gridcpu->ncorner_max[j] = gridcpu->ncorner_min[j] + gridcpu->nsize[j];
          gridcpu->gncorner_min[j] = gridcpu->ncorner_min[j] - Nghost[j] * grid->dn[j];
          gridcpu->gncorner_max[j] = gridcpu->ncorner_max[j] + Nghost[j] * grid->dn[j];
          gridcpu->gnsize[j] = gridcpu->gncorner_max[j] - gridcpu->gncorner_min[j];
          gridcpu->dn[j] = grid->dn[j];
          for (ii = INF; ii <= SUP; ii++)
          {
            if (grid->BoundaryConditions[j][ii] > 0)
              gridcpu->iface[j][ii] = grid->BoundaryConditions[j][ii];
            else
              gridcpu->iface[j][ii] = -1;
            /* -1: ghost in other tGrid, 0: ghost in same tGrid, >0: true BC */
            if (i[j] > 0)
              gridcpu->iface[j][INF] = 0; /* Means faces a brother (ie ghost in same tGrid) */
            if (i[j] < mx[j] - 1)
              gridcpu->iface[j][SUP] = 0;
          }
        }
        gridcpu->stride[0] = 1;
        gridcpu->stride[1] = (NDIM > 1 ? gridcpu->gncell[0] : 0);
        /* 0 is to avoid to call values outside of the mesh in the numerous
	   geometric specific terms in the hydro kernel */
        gridcpu->stride[2] = (NDIM > 2 ? gridcpu->gncell[0] * gridcpu->gncell[1] : 0);
        gridcpu->next = Grid_CPU_list;
        if (Grid_CPU_list != NULL)
          Grid_CPU_list->prev = gridcpu;
        gridcpu->prev = NULL;
        Grid_CPU_list = gridcpu;
        if (cpu == CPU_Rank)
          FillCPUGrid(gridcpu);
      }
    }
  }
}
