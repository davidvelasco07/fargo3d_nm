#define __NOPROTO
#define __GPU
#define THREADS 192
#include "fargo3d.h"


__global__ void kernel_Fluxes_X(real* flux, real* facefluxes_inf, real* facefluxes_sup, int size, int nx, int ny, int pitch, int stride, int FluxIndex) {
  int index, i, j, k;
  index=threadIdx.x+blockDim.x*blockIdx.x;
  
  if( index < size){
    k = index/ny;
    j = index - k*ny;
    i = NGHX + FLUXWIDTH;
    facefluxes_inf[j+k*ny+FluxIndex*size] += flux[i+(j+NGHY)*pitch+(k+NGHZ)*stride]; 
    i = NGHX + nx - FLUXWIDTH;
    facefluxes_sup[j+k*ny+FluxIndex*size] += flux[i+(j+NGHY)*pitch+(k+NGHZ)*stride]; 
  }
}
extern "C"
void Save_Face_Flux_X_gpu()
{
  tGrid_CPU *grid  = Flux->desc;
  int size = grid->ncell[_Y_]*grid->ncell[_Z_];
  kernel_Fluxes_X<<<(size+THREADS-1)/THREADS,THREADS>>>(Flux->field_gpu, Fluxes[_X_][0], Fluxes[_X_][1], size, grid->ncell[_X_], grid->ncell[_Y_], grid->Pitch_gpu, grid->Stride_gpu, FluxIndex);
  cudaDeviceSynchronize();
}

__global__ void kernel_Fluxes_Y(real* flux, real* facefluxes_inf, real* facefluxes_sup, int size, int nx, int ny, int pitch, int stride, int FluxIndex) {
  int index, i, j, k;
  index=threadIdx.x+blockDim.x*blockIdx.x;
  
  if( index < size){
    
    k = index/nx;
    i = index - k*nx;
    j = NGHY + FLUXWIDTH;
    facefluxes_inf[i+k*nx+FluxIndex*size] += flux[i+NGHX + j*pitch + (k+NGHZ)*stride]; 
    j = NGHY + ny - FLUXWIDTH;
    facefluxes_sup[i+k*nx+FluxIndex*size] += flux[i+NGHX + j*pitch + (k+NGHZ)*stride];
  
  }
}

extern "C"
void Save_Face_Flux_Y_gpu()
{
  tGrid_CPU *grid  = Flux->desc;
  int size = grid->ncell[_X_]*grid->ncell[_Z_];
  kernel_Fluxes_Y<<<(size+THREADS-1)/THREADS,THREADS>>>(Flux->field_gpu, Fluxes[_Y_][0], Fluxes[_Y_][1], size, grid->ncell[_X_], grid->ncell[_Y_], grid->Pitch_gpu, grid->Stride_gpu, FluxIndex);
  cudaDeviceSynchronize();
}

__global__ void kernel_Fluxes_Z(real* flux, real* facefluxes_inf, real* facefluxes_sup, int size, int nx, int nz, int pitch, int stride, int FluxIndex) {
  int index, i, j, k;
  index=threadIdx.x+blockDim.x*blockIdx.x;
  
  if( index < size){
    j = index/nx;
    i = index - j*nx;
    k = NGHZ + FLUXWIDTH;
    facefluxes_inf[i+j*nx+FluxIndex*size] += flux[i+NGHX + (j+NGHY)*pitch + k*stride]; 
    k = NGHZ + nz - FLUXWIDTH;
    facefluxes_sup[i+j*nx+FluxIndex*size] += flux[i+NGHX + (j+NGHY)*pitch + k*stride];
  }
}

extern "C"
void Save_Face_Flux_Z_gpu()
{
  tGrid_CPU *grid  = Flux->desc;
  int size = grid->ncell[_X_]*grid->ncell[_Y_];
  kernel_Fluxes_Z<<<(size+THREADS-1)/THREADS,THREADS>>>(Flux->field_gpu, Fluxes[_Z_][0], Fluxes[_Z_][1], size, grid->ncell[_X_], grid->ncell[_Z_], grid->Pitch_gpu, grid->Stride_gpu, FluxIndex);
  cudaDeviceSynchronize();
}




