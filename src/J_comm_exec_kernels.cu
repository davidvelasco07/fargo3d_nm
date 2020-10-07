#define __NOPROTO
#define __GPU
#include "fargo3d.h"

#define THREADS 128

__global__ void kernel_UP (real **fields, real *Ymin, real *Zmin, real OF, struct parms srcParms, real **buffers, struct gpucomm *CommArray, int *centered, int nfields, int number){

  int bufferid = blockIdx.x / nfields;
  struct gpucomm com = CommArray[bufferid];
  real *buffer = buffers[bufferid];
  int le = blockIdx.x - nfields*bufferid;
  int index, pitch, pitchd, stride, strided, size, sized, faceside, facedim;
  int id[3], imin[3], displacement[3], h;
  real s, vp, ym, sr;
#ifndef FASTCOMMGPU
  if ( number == com.srcg ){
#endif
    faceside = com.faceside;
    facedim =  com.facedim;
    size = com.size;
    sized = com.sizeD/sizeof(real);
    pitch = com.Pitch;
    stride = com.Stride;
    pitchd = com.pitchD/sizeof(real);
    strided = com.strideD/sizeof(real);
    imin[0] = com.xmin;
    imin[1] = com.ymin;
    imin[2] = com.zmin;
    
    
    for(h=0; h<3; h++)
      displacement[h]=0;
    
#if COMM == INTERNAL	
    displacement[facedim] = (1-faceside)*(1-centered[le*3+facedim]);
#endif
#if COMM == EXTERNAL
    displacement[facedim] = faceside*(1-centered[le*3+facedim]);
#endif
    
    for(index = threadIdx.x; index < size; index += blockDim.x){
      id[2] = index/stride;
      id[1] = (index - id[2]*stride)/pitch;
      id[0] = index - id[2]*stride - id[1]*pitch;    
      
      s = buffer[id[0] + id[1]*pitchd + id[2]*strided + le*sized];
      for (h = 0; h < 3; h++)
	id[h] += imin[h]+displacement[h];
#ifdef  COMMADAPT
      if(centered[le*3] == 0){
	ym = 0.5 *(Ymin[id[1]] + Ymin[id[1]+1]);
	vp = ym * sin( 0.5*(Zmin[id[2]]+Zmin[id[2]+1]))*OF;
	sr = sqrt(ym);
	s = s/sr - vp;
      }
#endif
      fields[le][id[0]+id[1]*srcParms.pitch+id[2]*srcParms.stride] = s;
    }
#ifndef FASTCOMMGPU
  }
#endif    
}

__global__ void kernel_DOWN (real **fields, struct parms dstParms, real **buffers, struct gpucommDown *CommArray, int *centered, int nfields, int number){

  int bufferid = blockIdx.x / nfields;
  struct gpucommDown com = CommArray[bufferid];
  real *buffer = buffers[bufferid];
  int le = blockIdx.x - nfields*bufferid;
  int index, pitch, pitchd, stride, strided, size, sized, faceside, facedim;
  int md, ms, id[3], imin[3], iface[3], displacement[3], cut[3], allcentered, stagdim, h;
  real s;
#ifndef FASTCOMMGPU
  if ( number == com.srcg ){
#endif
    faceside = com.faceside;
    facedim =  com.facedim;
    size = com.size;
    sized = com.sizeD/sizeof(real);
    pitch = com.Pitch;
    stride = com.Stride;
    pitchd = com.pitchD/sizeof(real);
    strided = com.strideD/sizeof(real);
    imin[0] = com.xmin;
    imin[1] = com.ymin;
    imin[2] = com.zmin;
    iface[0] = com.srcxface;
    iface[1] = com.srcyface;
    iface[2] = com.srczface;
    for(h=0; h<3; h++){
      displacement[h]=0;
      cut[h]=0;
    }
    displacement[facedim] = (1-faceside)*(1-centered[le*3+facedim]);
    allcentered = 1;
    stagdim=0;
    for(h=0; h<3; h++){
      allcentered *= centered[le*3+h];
      stagdim += (1-centered[le*3+h])*h;
    }
    if (allcentered == 0 && facedim != stagdim && iface[stagdim] == -1)cut[stagdim]=1;
    
    for(index = threadIdx.x; index < size; index += blockDim.x){
      id[2] = index/stride;
      id[1] = (index - id[2]*stride)/pitch;
      id[0] = index - id[2]*stride - id[1]*pitch;    
      ms = id[0] + id[1]*pitchd + id[2]*strided;
      if(id[stagdim] >= cut[stagdim]){
	s = buffer[ms + le*sized];
	for (h = 0; h < 3; h++)
	  id[h] += imin[h]+displacement[h];
	md = id[0]+id[1]*dstParms.pitch+id[2]*dstParms.stride;
	fields[le][md] = s;
	
      }
    }
#ifndef FASTCOMMGPU
  }
#endif
}


__global__ void kernel_RESETFLUX (real **facefluxes,  struct gpucommFlux *CommArray, int nfields){

  int bufferid = blockIdx.x / nfields;
  struct gpucommFlux com = CommArray[bufferid];
  int le = blockIdx.x - nfields*bufferid;
  int index, size, faceside, dim, dim1, n1, md;
  int id[3], Refine[3];
  faceside = com.faceside;
  size = com.sizeS;
  dim =  com.facedim;
  dim1 = (dim == 0);
  Refine[0] = com.xrefine+1;
  Refine[1] = com.yrefine+1;
  Refine[2] = com.zrefine+1;
  n1 = com.n1*Refine[dim1];

  
    for(index = threadIdx.x; index < size; index += blockDim.x){
  	id[1] = index/n1;
	id[0] = index - id[1]*n1;
      	md = id[0]+id[1]*n1;
	facefluxes[dim*2+faceside][md+le*size] = 0;
     }
 
    
}

extern "C"
void EXECUP (tGrid_CPU *grid, int fields)
{
  int number = -1;
#ifndef FASTCOMMGPU 
  jCommunicator *com;
  com = ComListGhost;
  while (com != NULL) {
    if ((com->destg == grid) && (grid->level > com->srcg->level) && (number != com->srcg->parent)) {
      number = com->srcg->parent;
      kernel_UP<<<grid->dst.nbCommsUp*fields,THREADS>>>(grid->dest, grid->Ymin_d, grid->Zmin_d, OMEGAFRAME, grid->gpuparms, grid->dst.BuffersUp, grid->dst.ParmsUp, grid->centered_gpu, fields, number);
      //cudaDeviceSynchronize();
    }
    com=com->next;
  }
#else
  if(grid->dst.nbCommsUp !=0){
    kernel_UP<<<grid->dst.nbCommsUp*fields,THREADS>>>(grid->dest, grid->Ymin_d, grid->Zmin_d, OMEGAFRAME, grid->gpuparms, grid->dst.BuffersUp, grid->dst.ParmsUp, grid->centered_gpu, fields, number);
    //cudaDeviceSynchronize();
  }
#endif
}		

extern "C"
void EXECDOWN (tGrid_CPU *grid, int fields)
{
  int number = -1;
#ifndef FASTCOMMGPU  
  jCommunicator *com;
  com = ComListMean;
  while (com != NULL) {
    if ((com->destg == grid) && (number != com->srcg->parent)) {
      number = com->srcg->parent;
      kernel_DOWN<<<grid->dst.nbCommsDown*fields,THREADS>>>(grid->dest,  grid->gpuparms, grid->dst.BuffersDown, grid->dst.ParmsDown, grid->centered_gpu, fields, number);
      //cudaDeviceSynchronize();
    }
    com=com->next;
  }
#else
  if(grid->dst.nbCommsDown !=0){
    kernel_DOWN<<<grid->dst.nbCommsDown*fields,THREADS>>>(grid->dest,  grid->gpuparms, grid->dst.BuffersDown, grid->dst.ParmsDown, grid->centered_gpu, fields, number);
    //cudaDeviceSynchronize();
  }
#endif
}		

extern "C"
void RESETFLUX (tGrid_CPU *grid, int fields)
{
  FluidPatch *fluid;
  fluid = grid->fluid;
  if(grid->src.nbCommsFlux !=0){
    while (fluid != NULL){
      kernel_RESETFLUX<<<grid->src.nbCommsFlux*fields,THREADS>>>(fluid->FluxesGPU, grid->src.ParmsFlux, fields);
      fluid = fluid->next;
    }
  }
}
