#define __NOPROTO
#define __GPU

#include "fargo3d.h"
#include <cudpp.h>

#define TILE 16

static CUDPPHandle theCudpp;
static CUDPPHandle scanplan = 0;
static CUDPPResult res;
static CUDPPConfiguration config;
static boolean PlanInitialized = NO;
static real *Transp_in;
static real *Transp_out;
static size_t tPitch;

void InitPlan () {
  // We need two arrays as scans on the GPU cannot be performed in place.
  cudaMallocPitch ((void **)&Transp_in,  &tPitch, (size_t)((Ny+2*NGHY)*sizeof(real)), (size_t)((Nx+2*NGHX)*(Nz+2*NGHZ)));
  check_errors ("First memory allocation on device in RT_IRR_LocalScan.cu");
  cudaMallocPitch ((void **)&Transp_out, &tPitch, (size_t)((Ny+2*NGHY)*sizeof(real)), (size_t)((Nx+2*NGHX)*(Nz+2*NGHZ)));
  check_errors ("Second memory allocation on device in RT_IRR_LocalScan.cu");
  tPitch /= sizeof(real);
  cudppCreate (&theCudpp);
  config.op = CUDPP_ADD;
#ifndef FLOAT
  config.datatype = CUDPP_DOUBLE;
#else
  config.datatype = CUDPP_FLOAT;
#endif
  config.algorithm = CUDPP_SCAN;
  config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_INCLUSIVE;
  res = cudppPlan (theCudpp, &scanplan, config, (int)(Ny+2*NGHY), (int)((Nx+2*NGHX)*(Nz+2*NGHZ)), tPitch);  
  if (res != CUDPP_SUCCESS) {
    fprintf (stderr, "Error creating CUDPPPlan in %s at line %d\n", __FILE__, __LINE__-2);
    exit (EXIT_FAILURE);
  }
  PlanInitialized = YES;
}
  
template <int direction>
__global__ void TransposeXY3D (real *in, real *out, int pitch_in, int pitch_out,\
			       int nx, int ny, int nz) 
{ //ny = Ny+2*NGHY in direct transpose.
  __shared__ real localtile[TILE*(TILE+1)]; //+1 to avoid bank conflicts
  int xg = blockIdx.x * TILE + threadIdx.x;
  int yg = blockIdx.y * TILE + threadIdx.y;
  int zg = blockIdx.z;
  int index;
  if ((xg < nx) && (yg < ny)) {
    index = xg + (yg + ny*zg)*pitch_in;
    //The test below is non divergent
    if ((yg < NGHY) && (direction == 1)) // In the first transpose
					 // (prior to the scan) it is
					 // important to have zeroes
					 // in the inner ghosts in
					 // radius.
      localtile[threadIdx.x+(TILE+1)*threadIdx.y] = 0.0;
    else
      localtile[threadIdx.x+(TILE+1)*threadIdx.y] = in[index];
  }
  __syncthreads();
  //We calculate the new global indices
  xg = blockIdx.y * TILE + threadIdx.x; //ensures coalesced write
  yg = blockIdx.x * TILE + threadIdx.y;
  if ((xg < ny) && (yg < nx)) {
    index = xg + (yg + nx*zg)*pitch_out;
    out[index] = localtile[threadIdx.y+(TILE+1)*threadIdx.x]; // Bank conflicts free
  }
}

extern "C" void InProcessRadialScan_gpu (Field *Array) {
  if (PlanInitialized == NO) InitPlan();
  
  INPUT (Array);
  OUTPUT (Array);
  dim3 grid, block;
  
  block.x = block.y = TILE;
  block.z = 1;
  grid.x = (Nx+2*NGHX+TILE-1)/TILE;
  grid.y = (Ny+2*NGHY+TILE-1)/TILE;
  grid.z = Nz+2*NGHZ;
  TransposeXY3D<+1><<<grid, block>>>(Array->field_gpu, Transp_in, Pitch_gpu, tPitch, Nx+2*NGHX, Ny+2*NGHY, Nz+2*NGHZ);
  check_errors ("First transpose in RT_IRR_LocalScan.cu");

  cudppMultiScan (scanplan, Transp_out, Transp_in, Ny+2*NGHY, (Nx+2*NGHX)*(Nz+2*NGHZ));

  grid.x = (Ny+2*NGHY+TILE-1)/TILE;
  grid.y = (Nx+2*NGHX+TILE-1)/TILE;
  TransposeXY3D<-1><<<grid, block>>>(Transp_out, Array->field_gpu, tPitch, Pitch_gpu, Ny+2*NGHY, Nx+2*NGHX, Nz+2*NGHZ);
  check_errors ("Second transpose in RT_IRR_LocalScan.cu");
}
