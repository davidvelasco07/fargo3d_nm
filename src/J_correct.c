#include "fargo3d.h"

void OverWriteBoundaryFluxes_cpu(int dim)
{
  //<USER_DEFINED>
  OUTPUT(Flux);
  //<\USER_DEFINED>
  long i, side, dim1, dim2;
  long imin[3], imax[3], size[3], j, k, m, le, n;
  jCommunicator *com;
  real *dest[8][3], *buffer;
  real *flux = Flux->field_cpu;
  com = ComListFlux;

  while (com != NULL)
  {
    if ((com->destg == Current_Grid) && (dim == com->facedim))
    {
      side = com->faceside;
      dim1 = (dim == 0);
      dim2 = 2 - (dim == 2);
      for (i = 0; i < 3; i++)
      {
        imin[i] = com->imin_dest[i];
        imax[i] = com->imax_dest[i];
        size[i] = imax[i] - imin[i];
      }

      buffer = com->buffer;
      i = (side == INF ? imax[dim] : imin[dim]);
      for (j = imin[dim1]; j < imax[dim1]; j++)
      {
        for (k = imin[dim2]; k < imax[dim2]; k++)
        {
          m = i * jstride[dim] + j * jstride[dim1] + k * jstride[dim2];
          n = (j - imin[dim1]) + (k - imin[dim2]) * size[dim1];
          if (FirstPassXAdvection)
            buffer[n + FluxIndex * size[dim1] * size[dim2]] -= flux[m];
          else
            flux[m] = buffer[n + FluxIndex * size[dim1] * size[dim2]];
        }
      }
    }

    com = com->next;
  }
  /* The FluxIndex variable allows to store the fluxes of different
     quantities in the same array. It should only be incremented in
     this function. Both the sender and receiver (fine and coarse
     level) go through the line below, regardless of whether they have
     a finer level or not. This line increments the offset of the
     quantity every time a new quantity is transported. FluxIndex
     should be reset to 0 every time we change dimension. */
  FluxIndex++;
}

void OverWriteBoundaryFluxes(int dim){
#ifndef GPU
	OverWriteBoundaryFluxes_cpu(dim);
#else
#ifndef COMMGPU
	OverWriteBoundaryFluxes_cpu(dim);
#else
	Output_GPU(Flux, __LINE__, __FILE__);
  OverWriteBoundaryFluxes_gpu(dim);
#endif
#endif
}