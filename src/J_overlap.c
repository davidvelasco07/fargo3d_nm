#include "fargo3d.h"

void SetOverlapFlag()
{
  long size, i, j, k, diff, m;
  long rcor[6], gcor[6], inter[6];
  long imin[3], imax[3];
  long gncell[3], stride[3];
  tGrid_CPU *g, *r;
  char *flags;
  int nx, ny, nz, N;
  g = Grid_CPU_list;
  while (g != NULL)
  {
    getgridsize(g, gncell, stride);
    if (g->cpu == CPU_Rank)
    {
      size = 1;
      for (i = 0; i < NDIM; i++)
        size *= gncell[i];
      flags = (char *)prs_malloc(size * sizeof(char));
      for (i = 0; i < size; i++)
        flags[i] = 0;
      g->Hidden = flags;
      r = Grid_CPU_list;
      diff = 1;
      while (r != NULL)
      {
        if (r == g->next)
          diff = 0; /* Regardless of their respective
				     levels, g is in the forefront of
				     subsequent grids */
        if (r->level - g->level >= diff)
        {
          for (i = 0; i < 3; i++)
          {
            rcor[i] = r->ncorner_min[i];
            rcor[i + 3] = r->ncorner_max[i];
            gcor[i] = g->ncorner_min[i];
            gcor[i + 3] = g->ncorner_max[i];
          }
          if (CubeIntersect(rcor, gcor, inter))
          {
            //printf("grid %d of level %d is overlaped by grid %d of level %d\n",g->number, g->level, r->number, r->level);
            for (i = 0; i < 3; i++)
            {
              imin[i] = (inter[i] - gcor[i]) / (g->Parent->dn[i]) + Nghost[i];
              imax[i] = (inter[i + 3] - gcor[i]) / (g->Parent->dn[i]) + Nghost[i];
            }
            for (i = imin[0]; i < imax[0]; i++)
            {
              for (j = imin[1]; j < imax[1]; j++)
              {
                for (k = imin[2]; k < imax[2]; k++)
                {
                  m = i * stride[0] + j * stride[1] + k * stride[2];
                  flags[m] = 1;
                }
              }
            }
          }
        }
        r = r->next;
      }
#ifdef GPU
      masterprint("Allocating and copying Hidden to device\n");
      cudaMalloc((void **)&g->Hidden_d, Maxsize_gpu * sizeof(char));
      nx = gncell[0];
      ny = gncell[1];
      nz = gncell[2];
      if (nx == 1)
        cudaMemcpy2D(g->Hidden_d, g->Stride_gpu * sizeof(char), g->Hidden, ny * sizeof(char), ny * sizeof(char), nz, cudaMemcpyHostToDevice);
      else
        cudaMemcpy2D(g->Hidden_d, g->Pitch_gpu * sizeof(char), g->Hidden, nx * sizeof(char), nx * sizeof(char), ny * nz, cudaMemcpyHostToDevice);
      check_errors("cudaMemcpy2D");
#endif
    }
    g = g->next;
  }
}
