#include "fargo3d.h"

void ResetFluxesLevel(int lev)
{
  int dim, side, dimp1, dimp2, k, size, nvar;
  tGrid_CPU *item;
  FluidPatch *fluid;
  item = Grid_CPU_list;
  nvar = 1 + 2 * NDIM;
#ifdef ADIABATIC
  nvar++;
#endif
  do
  {
    if (item->cpu == CPU_Rank)
    {
      if (item->level == lev)
      {
#ifndef GPUCOMM
        fluid = item->fluid;
        do{
          for (dim = 0; dim < 3; dim++)
          {
            for (side = INF; side <= SUP; side++)
            {
              dimp1 = (dim == 0);
              dimp2 = 2 - (dim == 2);
              size = (item->ncell[dimp1]) * (item->ncell[dimp2]) * nvar;
              for (k = 0; k < size; k++)
                fluid->Fluxes[dim][side][k] = 0.0;
            }
          }
          fluid = fluid->next;
        }while(fluid != NULL);

#else
        RESETFLUX(item, nvar);
#endif
      }
    }
    item = item->next;
  } while (item != NULL);
}
