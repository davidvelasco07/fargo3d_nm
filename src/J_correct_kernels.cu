#define __NOPROTO
#define __GPU

#include "fargo3d.h"

#define THREADS 192

__global__ void kernel_OWFLUX (real *flux, struct parms srcParms, real **buffers, struct gpucommFlux *CommArray, int fluxindex, int dim, int number, boolean FirstPassXAdvection){

  struct gpucommFlux com = CommArray[blockIdx.x];
  real *buffer = buffers[blockIdx.x];
  int index, pitchd, strided, size, sized, facedim, faceside, n1, md, ms;
  int id[3], is[3], imin[3], h, dim1, dim2;
  real s;
  facedim =  com.facedim;
  faceside =  com.faceside;
  size = com.size;
  sized = com.sizeD/sizeof(real);
  n1 = com.n1;
  pitchd = com.pitchD/sizeof(real);
  strided = com.strideD/sizeof(real);
  dim1 = (facedim == 0);
  dim2 = 2 - (facedim == 2);
  imin[0] = com.xmin;
  imin[1] = com.ymin;
  imin[2] = com.zmin;
#ifndef FASTCOMMGPU
  if (com.srcg == number){
#endif
  if(facedim == dim){
    for(index = threadIdx.x; index < size; index += blockDim.x){
      is[dim2] = index/n1;
      is[dim1] = index - is[dim2]*n1;
      is[dim] = 0;
      ms = is[0]+is[1]*pitchd+is[2]*strided;
      if (FirstPassXAdvection) {
	for (h = 0; h < 3; h++)
	  id[h] = imin[h]+is[h];
	id[dim] += 1-faceside;
	md = id[0] + id[1]*srcParms.pitch + id[2]*srcParms.stride;
	buffer[ms+sized*fluxindex] -= flux[md];
      } else {
	s = buffer[ms+sized*fluxindex];
	for (h = 0; h < 3; h++)
	  id[h] = imin[h]+is[h];
	id[dim] += 1-faceside;
	md = id[0] + id[1]*srcParms.pitch + id[2]*srcParms.stride;
	flux[md]=s;
      }
    }
  }
#ifndef FASTCOMMGPU
  }
#endif   
}

extern "C"
void OverWriteBoundaryFluxes_gpu (int dim)
{
  tGrid_CPU *grid = Current_Grid;
  int number = -1;
#ifndef FASTCOMMGPU 
  jCommunicator *com;
  com = ComListFlux;
  while (com != NULL) {
    if ((com->destg == Current_Grid) && (dim == com->facedim) && (number != com->srcg->parent)) {
      number = com->srcg->parent;
      kernel_OWFLUX<<<grid->dst.nbCommsFlux,THREADS>>>(Flux->field_gpu, grid->gpuparms, grid->dst.BuffersFlux, grid->dst.ParmsFlux, FluxIndex, dim, number, FirstPassXAdvection);
      cudaThreadSynchronize();
    }
    com=com->next;
  }
#else
  if(grid->dst.nbCommsFlux !=0){ 
    kernel_OWFLUX<<<grid->dst.nbCommsFlux,THREADS>>>(Flux->field_gpu, grid->gpuparms, grid->dst.BuffersFlux, grid->dst.ParmsFlux, FluxIndex, dim, number, FirstPassXAdvection);
    cudaThreadSynchronize();
  }
#endif
  FluxIndex++;
}

