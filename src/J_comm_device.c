#include "fargo3d.h"

extern MPI_Request *Req;
#ifdef GPU
struct cudaPitchedPtr source_gpu[80];
struct cudaPitchedPtr dest_gpu[80];
#endif
void JCommunicatorGPU(jCommunicator *comm)
{
#ifdef GPU
  size_t pitch;
  struct cudaPitchedPtr buffer;
  struct cudaExtent extent;
  int nvar = 2 + NDIM;
  if (comm->type == FLUX)
  {
    nvar = 1 + 2 * NDIM;
#ifdef ADIABATIC
    nvar++;
#endif
  }
  int xminsrc = comm->imin_src[0];
  int yminsrc = comm->imin_src[1];
  int zminsrc = comm->imin_src[2];
  int xmaxsrc = comm->imax_src[0];
  int ymaxsrc = comm->imax_src[1];
  int zmaxsrc = comm->imax_src[2];
  int xmindest = comm->imin_dest[0];
  int ymindest = comm->imin_dest[1];
  int zmindest = comm->imin_dest[2];
  int xmaxdest = comm->imax_dest[0];
  int ymaxdest = comm->imax_dest[1];
  int zmaxdest = comm->imax_dest[2];
  comm->OnSrc.srcArray = NULL;
  comm->OnSrc.dstArray = NULL;
  comm->OnDst.srcArray = NULL;
  comm->OnDst.dstArray = NULL;
  comm->OnSrc.srcPos = make_cudaPos(xminsrc * sizeof(real), yminsrc, zminsrc);
  comm->OnDst.dstPos = make_cudaPos(xmindest * sizeof(real), ymindest, zmindest);
  comm->OnSrc.extent = make_cudaExtent((xmaxdest - xmindest) * sizeof(real), ymaxdest - ymindest, zmaxdest - zmindest);
  comm->OnDst.extent = make_cudaExtent((xmaxdest - xmindest) * sizeof(real), ymaxdest - ymindest, zmaxdest - zmindest);
  comm->OnSrc.kind = cudaMemcpyDeviceToDevice;
  comm->OnDst.kind = cudaMemcpyDeviceToDevice;
  extent = make_cudaExtent((xmaxdest - xmindest) * sizeof(real), (ymaxdest - ymindest), (zmaxdest - zmindest) * nvar);
  cudaMalloc3D(&buffer, extent);
  cudaMemset3D(buffer, 1, extent);
  check_errors("allocating comm buffer on GPU");
  comm->OnSrc.dstPtr = buffer;
  comm->OnDst.srcPtr = buffer;
  comm->bufferGPU = buffer;
  comm->dz = zmaxdest - zmindest;
  comm->yzsize = (ymaxdest - ymindest) * (zmaxdest - zmindest);
#endif
}

void CPU_Grid_InitChainedLists()
{
#ifdef GPU
  tGrid_CPU *item = Grid_CPU_list;
  while (item != NULL)
  {
    item->src.CommListDown = NULL;
    item->src.CommListUp = NULL;
    item->src.CommListSame = NULL;
    item->src.CommListFlux = NULL;
    item->dst.CommListDown = NULL;
    item->dst.CommListUp = NULL;
    item->dst.CommListSame = NULL;
    item->dst.CommListFlux = NULL;
    item = item->next;
  }
#endif
}

void Comm_SubList()
{
#ifdef GPU
  jCommunicator *com;
  CommHash *hash;
  CPU_Grid_InitChainedLists();
  com = ComListMean;
  while (com != NULL)
  {
    if (com->CPU_src == CPU_Rank)
    {
      hash = (CommHash *)prs_malloc(sizeof(CommHash));
      hash->next = com->srcg->src.CommListDown;
      hash->com = com;
      com->srcg->src.CommListDown = hash;
      com->srcg->src.nbCommsDown += 1;
    }
    if (com->CPU_dest == CPU_Rank)
    {
      hash = (CommHash *)prs_malloc(sizeof(CommHash));
      hash->next = com->destg->dst.CommListDown;
      hash->com = com;
      com->destg->dst.CommListDown = hash;
      com->destg->dst.nbCommsDown += 1;
    }
    com = com->next;
  }

  com = ComListGhost;
  while (com != NULL)
  {
    if (com->CPU_src == CPU_Rank)
    {
      if (com->dest_level > com->src_level)
      {
        hash = (CommHash *)prs_malloc(sizeof(CommHash));
        hash->next = com->srcg->src.CommListUp;
        hash->com = com;
        com->srcg->src.CommListUp = hash;
        com->srcg->src.nbCommsUp += 1;
      }
      else if (com->dest_level == com->src_level)
      {
        hash = (CommHash *)prs_malloc(sizeof(CommHash));
        hash->next = com->srcg->src.CommListSame;
        hash->com = com;
        com->srcg->src.CommListSame = hash;
      }
    }
    if (com->CPU_dest == CPU_Rank)
    {
      if (com->dest_level > com->src_level)
      {
        hash = (CommHash *)prs_malloc(sizeof(CommHash));
        hash->next = com->destg->dst.CommListUp;
        hash->com = com;
        com->destg->dst.CommListUp = hash;
        com->destg->dst.nbCommsUp += 1;
      }
      else if (com->dest_level == com->src_level)
      {
        hash = (CommHash *)prs_malloc(sizeof(CommHash));
        hash->next = com->destg->dst.CommListSame;
        hash->com = com;
        com->destg->dst.CommListSame = hash;
      }
    }
    com = com->next;
  }
  com = ComListFlux;
  while (com != NULL)
  {
    if (com->CPU_src == CPU_Rank)
    {
      hash = (CommHash *)prs_malloc(sizeof(CommHash));
      hash->next = com->srcg->src.CommListFlux;
      hash->com = com;
      com->srcg->src.CommListFlux = hash;
      com->srcg->src.nbCommsFlux += 1;
    }
    if (com->CPU_dest == CPU_Rank)
    {
      hash = (CommHash *)prs_malloc(sizeof(CommHash));
      hash->next = com->destg->dst.CommListFlux;
      hash->com = com;
      com->destg->dst.CommListFlux = hash;
      com->destg->dst.nbCommsFlux += 1;
    }
    com = com->next;
  }
  FillCommArrays();
  check_errors("FillCommArrays");
#endif
}

void FillCommArrays()
{
#ifdef GPU
  tGrid_CPU *grid;
  grid = Grid_CPU_list;
  CommHash *hash;
  int i, j, h, d, le, dimp1, dimp2;
  void *ptr_gpu;
  struct gpucomm ArrayUp[100];
  struct gpucommDown ArrayDown[100];
  struct gpucommFlux ArrayFlux[100];
  real *source[100];
  real *buffers[100];
  int field, nvar;
  int centered[60];
  FluidPatch *fluid;
  int fieldtype[20];
  int options = DENS | ENERGY | VX | VY | VZ;

  nvar = BuildFieldType(fieldtype, options);
  while (grid != NULL)
  {
    if (grid->cpu == CPU_Rank)
    {
#ifdef X
      grid->gpuparms.pitch = grid->Pitch_gpu;
      grid->gpuparms.stride = grid->Stride_gpu;
#else
      grid->gpuparms.pitch = grid->Stride_gpu;
      grid->gpuparms.stride = grid->Stride_gpu * grid->gncell[1];
#endif
      grid->gpuparms.cxmin = grid->gncorner_min[0];
      grid->gpuparms.cymin = grid->gncorner_min[1];
      grid->gpuparms.czmin = grid->gncorner_min[2];
      grid->gpuparms.dx = grid->dn[0];
      grid->gpuparms.dy = grid->dn[1];
      grid->gpuparms.dz = grid->dn[2];
      grid->gpuparms.sizex = grid->ncell[0];
      grid->gpuparms.sizey = grid->ncell[1];
      grid->gpuparms.sizez = grid->ncell[2];

      cudaMalloc((void **)&grid->centered_gpu, nvar * sizeof(int));
      cudaMalloc((void **)&grid->source, nvar * sizeof(real));
      cudaMalloc((void **)&grid->dest, nvar * sizeof(real));

      hash = grid->src.CommListUp;
      i = 0;
      while (hash != NULL)
      {
        ArrayUp[i].facedim = hash->com->facedim;
        ArrayUp[i].faceside = hash->com->faceside;
        ArrayUp[i].size = hash->com->OnSrc.extent.width * hash->com->OnSrc.extent.height * hash->com->OnSrc.extent.depth / sizeof(real);
        ArrayUp[i].pitchD = hash->com->OnSrc.dstPtr.pitch;
        ArrayUp[i].strideD = hash->com->OnSrc.extent.height * hash->com->OnSrc.dstPtr.pitch;
        ArrayUp[i].sizeD = hash->com->OnSrc.extent.height * hash->com->OnSrc.dstPtr.pitch * hash->com->dz;
        ArrayUp[i].Pitch = hash->com->OnSrc.extent.width / sizeof(real);
        ArrayUp[i].Stride = hash->com->OnSrc.extent.width * hash->com->OnSrc.extent.height / sizeof(real);
        ArrayUp[i].xmin = hash->com->imin_dest[0];
        ArrayUp[i].ymin = hash->com->imin_dest[1];
        ArrayUp[i].zmin = hash->com->imin_dest[2];
        ArrayUp[i].cxmin = hash->com->destg->gncorner_min[0];
        ArrayUp[i].cymin = hash->com->destg->gncorner_min[1];
        ArrayUp[i].czmin = hash->com->destg->gncorner_min[2];
        ArrayUp[i].dx = (int)hash->com->destg->dn[0];
        ArrayUp[i].dy = (int)hash->com->destg->dn[1];
        ArrayUp[i].dz = (int)hash->com->destg->dn[2];
        buffers[i] = hash->com->bufferGPU.ptr;
        i++;
        hash = hash->next;
      }
      cudaMalloc((void **)&ptr_gpu, (grid->src.nbCommsUp) * sizeof(struct gpucomm));
      grid->src.ParmsUp = (struct gpucomm *)ptr_gpu;
      cudaMemcpy(grid->src.ParmsUp, &ArrayUp, (grid->src.nbCommsUp) * sizeof(struct gpucomm), cudaMemcpyHostToDevice);
      cudaMalloc((void **)&grid->src.BuffersUp, (grid->src.nbCommsUp) * sizeof(real));
      cudaMemcpy(grid->src.BuffersUp, &buffers, (grid->src.nbCommsUp) * sizeof(real), cudaMemcpyHostToDevice);

      hash = grid->src.CommListDown;
      i = 0;
      while (hash != NULL)
      {
        ArrayDown[i].facedim = hash->com->facedim;
        ArrayDown[i].faceside = hash->com->faceside;
        ArrayDown[i].size = hash->com->OnSrc.extent.width * hash->com->OnSrc.extent.height * hash->com->OnSrc.extent.depth / sizeof(real);
        ArrayDown[i].pitchD = hash->com->OnSrc.dstPtr.pitch;
        ArrayDown[i].strideD = hash->com->OnSrc.extent.height * hash->com->OnSrc.dstPtr.pitch;
        ArrayDown[i].sizeD = hash->com->OnSrc.extent.height * hash->com->OnSrc.dstPtr.pitch * hash->com->dz;
        ArrayDown[i].Pitch = hash->com->OnSrc.extent.width / sizeof(real);
        ArrayDown[i].Stride = hash->com->OnSrc.extent.width * hash->com->OnSrc.extent.height / sizeof(real);
        ArrayDown[i].xmin = hash->com->imin_src[0];
        ArrayDown[i].ymin = hash->com->imin_src[1];
        ArrayDown[i].zmin = hash->com->imin_src[2];
        ArrayDown[i].xrefine = Refine[0];
        ArrayDown[i].yrefine = Refine[1];
        ArrayDown[i].zrefine = Refine[2];
        buffers[i] = hash->com->bufferGPU.ptr;
        i++;
        hash = hash->next;
      }
      cudaMalloc((void **)&ptr_gpu, (grid->src.nbCommsDown) * sizeof(struct gpucommDown));
      grid->src.ParmsDown = (struct gpucommDown *)ptr_gpu;
      cudaMemcpy(grid->src.ParmsDown, &ArrayDown, (grid->src.nbCommsDown) * sizeof(struct gpucommDown), cudaMemcpyHostToDevice);
      cudaMalloc((void **)&grid->src.BuffersDown, (grid->src.nbCommsDown) * sizeof(real));
      cudaMemcpy(grid->src.BuffersDown, &buffers, (grid->src.nbCommsDown) * sizeof(real), cudaMemcpyHostToDevice);

      hash = grid->src.CommListFlux;
      i = 0;
      while (hash != NULL)
      {
        dimp1 = (hash->com->facedim == 0);
        dimp2 = 2 - (hash->com->facedim == 2);
        ArrayFlux[i].facedim = hash->com->facedim;
        ArrayFlux[i].faceside = hash->com->faceside;
        ArrayFlux[i].size = hash->com->OnSrc.extent.width * hash->com->OnSrc.extent.height * hash->com->OnSrc.extent.depth / sizeof(real);
        ArrayFlux[i].pitchD = hash->com->OnSrc.dstPtr.pitch;
        ArrayFlux[i].strideD = hash->com->OnSrc.extent.height * hash->com->OnSrc.dstPtr.pitch;
        ArrayFlux[i].sizeD = hash->com->OnSrc.extent.height * hash->com->OnSrc.dstPtr.pitch * hash->com->dz;
        ArrayFlux[i].n1 = hash->com->imax_dest[dimp1] - hash->com->imin_dest[dimp1];
        ArrayFlux[i].xmin = hash->com->imin_src[0] - Nghost[0];
        ArrayFlux[i].ymin = hash->com->imin_src[1] - Nghost[1];
        ArrayFlux[i].zmin = hash->com->imin_src[2] - Nghost[2];
        ArrayFlux[i].xrefine = Refine[0];
        ArrayFlux[i].yrefine = Refine[1];
        ArrayFlux[i].zrefine = Refine[2];
        ArrayFlux[i].sizeS = grid->gncell[dimp1] * grid->gncell[dimp2];
        buffers[i] = hash->com->bufferGPU.ptr;
        i++;
        hash = hash->next;
      }
      cudaMalloc((void **)&ptr_gpu, (grid->src.nbCommsFlux) * sizeof(struct gpucommDown));
      grid->src.ParmsFlux = (struct gpucommFlux *)ptr_gpu;
      cudaMemcpy(grid->src.ParmsFlux, &ArrayFlux, (grid->src.nbCommsFlux) * sizeof(struct gpucommFlux), cudaMemcpyHostToDevice);
      cudaMalloc((void **)&grid->src.BuffersFlux, (grid->src.nbCommsFlux) * sizeof(real));
      cudaMemcpy(grid->src.BuffersFlux, &buffers, (grid->src.nbCommsFlux) * sizeof(real), cudaMemcpyHostToDevice);

      hash = grid->dst.CommListUp;
      i = 0;
      while (hash != NULL)
      {
        ArrayUp[i].facedim = hash->com->facedim;
        ArrayUp[i].faceside = hash->com->faceside;
        ArrayUp[i].size = hash->com->OnDst.extent.width * hash->com->OnDst.extent.height * hash->com->OnDst.extent.depth / sizeof(real);
        ArrayUp[i].pitchD = hash->com->OnDst.srcPtr.pitch;
        ArrayUp[i].strideD = hash->com->OnDst.extent.height * hash->com->OnDst.srcPtr.pitch;
        ArrayUp[i].sizeD = hash->com->OnDst.extent.height * hash->com->OnDst.srcPtr.pitch * hash->com->dz;
        ArrayUp[i].Pitch = hash->com->OnDst.extent.width / sizeof(real);
        ArrayUp[i].Stride = hash->com->OnDst.extent.width * hash->com->OnDst.extent.height / sizeof(real);
        ArrayUp[i].xmin = hash->com->imin_dest[0];
        ArrayUp[i].ymin = hash->com->imin_dest[1];
        ArrayUp[i].zmin = hash->com->imin_dest[2];
        ArrayUp[i].srcg = hash->com->srcg->parent;
        buffers[i] = hash->com->bufferGPU.ptr;
        i++;
        hash = hash->next;
      }
      cudaMalloc((void **)&ptr_gpu, (grid->dst.nbCommsUp) * sizeof(struct gpucomm));
      grid->dst.ParmsUp = (struct gpucomm *)ptr_gpu;
      cudaMemcpy(grid->dst.ParmsUp, &ArrayUp, (grid->dst.nbCommsUp) * sizeof(struct gpucomm), cudaMemcpyHostToDevice);
      cudaMalloc((void **)&grid->dst.BuffersUp, (grid->dst.nbCommsUp) * sizeof(real));
      cudaMemcpy(grid->dst.BuffersUp, &buffers, (grid->dst.nbCommsUp) * sizeof(real), cudaMemcpyHostToDevice);

      hash = grid->dst.CommListDown;
      i = 0;
      while (hash != NULL)
      {
        ArrayDown[i].facedim = hash->com->facedim;
        ArrayDown[i].faceside = hash->com->faceside;
        ArrayDown[i].size = hash->com->OnDst.extent.width * hash->com->OnDst.extent.height * hash->com->OnDst.extent.depth / sizeof(real);
        ArrayDown[i].pitchD = hash->com->OnDst.srcPtr.pitch;
        ArrayDown[i].strideD = hash->com->OnDst.extent.height * hash->com->OnDst.srcPtr.pitch;
        ArrayDown[i].sizeD = hash->com->OnDst.extent.height * hash->com->OnDst.srcPtr.pitch * hash->com->dz;
        ArrayDown[i].Pitch = hash->com->OnDst.extent.width / sizeof(real);
        ArrayDown[i].Stride = hash->com->OnDst.extent.width * hash->com->OnDst.extent.height / sizeof(real);
        ArrayDown[i].xmin = hash->com->imin_dest[0];
        ArrayDown[i].ymin = hash->com->imin_dest[1];
        ArrayDown[i].zmin = hash->com->imin_dest[2];
        ArrayDown[i].srcxface = hash->com->srcg->iface[0][INF];
        ArrayDown[i].srcyface = hash->com->srcg->iface[1][INF];
        ArrayDown[i].srczface = hash->com->srcg->iface[2][INF];
        ArrayDown[i].srcg = hash->com->srcg->parent;
        buffers[i] = hash->com->bufferGPU.ptr;
        i++;
        hash = hash->next;
      }
      cudaMalloc((void **)&ptr_gpu, (grid->dst.nbCommsDown) * sizeof(struct gpucommDown));
      grid->dst.ParmsDown = (struct gpucommDown *)ptr_gpu;
      cudaMemcpy(grid->dst.ParmsDown, &ArrayDown, (grid->dst.nbCommsDown) * sizeof(struct gpucommDown), cudaMemcpyHostToDevice);
      cudaMalloc((void **)&grid->dst.BuffersDown, (grid->dst.nbCommsDown) * sizeof(real));
      cudaMemcpy(grid->dst.BuffersDown, &buffers, (grid->dst.nbCommsDown) * sizeof(real), cudaMemcpyHostToDevice);

      hash = grid->dst.CommListFlux;
      i = 0;
      while (hash != NULL)
      {
        dimp1 = (hash->com->facedim == 0);
        dimp2 = 2 - (hash->com->facedim == 2);
        ArrayFlux[i].facedim = hash->com->facedim;
        ArrayFlux[i].faceside = hash->com->faceside;
        ArrayFlux[i].size = hash->com->OnDst.extent.width * hash->com->OnDst.extent.height * hash->com->OnDst.extent.depth / sizeof(real);
        ArrayFlux[i].pitchD = hash->com->OnDst.srcPtr.pitch;
        ArrayFlux[i].strideD = hash->com->OnDst.extent.height * hash->com->OnDst.srcPtr.pitch;
        ArrayFlux[i].sizeD = hash->com->OnDst.extent.height * hash->com->OnDst.srcPtr.pitch * hash->com->dz;
        ArrayFlux[i].n1 = hash->com->imax_dest[dimp1] - hash->com->imin_dest[dimp1];
        ArrayFlux[i].xmin = hash->com->imin_dest[0];
        ArrayFlux[i].ymin = hash->com->imin_dest[1];
        ArrayFlux[i].zmin = hash->com->imin_dest[2];
        ArrayFlux[i].xrefine = Refine[0];
        ArrayFlux[i].yrefine = Refine[1];
        ArrayFlux[i].zrefine = Refine[2];
        ArrayFlux[i].srcg = hash->com->srcg->parent;
        buffers[i] = hash->com->bufferGPU.ptr;
        i++;
        hash = hash->next;
      }
      cudaMalloc((void **)&ptr_gpu, (grid->dst.nbCommsFlux) * sizeof(struct gpucommDown));
      grid->dst.ParmsFlux = (struct gpucommFlux *)ptr_gpu;
      cudaMemcpy(grid->dst.ParmsFlux, &ArrayFlux, (grid->dst.nbCommsFlux) * sizeof(struct gpucommFlux), cudaMemcpyHostToDevice);
      cudaMalloc((void **)&grid->dst.BuffersFlux, (grid->dst.nbCommsFlux) * sizeof(real));
      cudaMemcpy(grid->dst.BuffersFlux, &buffers, (grid->dst.nbCommsFlux) * sizeof(real), cudaMemcpyHostToDevice);
    }
    grid = grid->next;
  }
  check_errors("FillArrays");

#endif
}

void ExecCommSameVar_gpu(int lev, int nvar, int *fieldtype)
{
#ifdef GPU
  /* Execute intra-level communications, for level lev. These are
     GHOST type communications */
  MPI_Request reqs[8], reqr[8];
  jCommunicator *comm;
  long i, j, k, le, m, n, d;
  int field, nbreqs = 0, nbreqr = 0;
  FluidPatch *fluid;
  //struct cudaPitchedPtr source[80];
  comm = ComListGhost;
  while (comm != NULL)
  {
    if ((comm->dest_level == lev) && (comm->src_level == lev))
    {
      if (comm->CPU_src == CPU_Rank)
      {
        comm->OnSrc.dstPtr = comm->bufferGPU;
        fluid = comm->srcg->fluid;
        BuildCommSource_gpu(fluid, fieldtype, nvar, 0);
        for (le = 0; le < nvar; le++)
        {
          comm->OnSrc.dstPos = make_cudaPos(0, 0, le * comm->dz);
          comm->OnSrc.srcPtr = source_gpu[le];
          cudaMemcpy3D(&(comm->OnSrc));
          check_errors("cudaMemcpy3D");
        }
      }
    }
    comm = comm->next;
  }
  FARGO_SAFE(ExecComm_gpu(ComListGhost, lev, lev, nvar, fieldtype));
  comm = ComListGhost;
  while (comm != NULL)
  {
    if ((comm->dest_level == lev) && (comm->src_level == lev))
    {
      if (comm->CPU_dest == CPU_Rank)
      {
        comm->OnDst.srcPtr = comm->bufferGPU;
        fluid = comm->destg->fluid;
        BuildCommDest_gpu(fluid, fieldtype, nvar, 0);
        for (le = 0; le < nvar; le++)
        {
          comm->OnDst.srcPos = make_cudaPos(0, 0, le * comm->dz);
          comm->OnDst.dstPtr = dest_gpu[le];
          cudaMemcpy3D(&(comm->OnDst));
          check_errors("cudaMemcpy3D");
        }
      }
    }
    comm = comm->next;
  }
  //FARGO_SAFE(ExecCommS_gpu (ComListGhost, lev, lev, nvar, fieldtype));
  //check_errors ("ExecCommS_gpu");
#endif
}

void ExecCommUpVar_gpu(int level, int nvar, int *fieldtype)
{
#ifdef GPU
  tGrid_CPU *grid;
  grid = Grid_CPU_list;
  FluidPatch *fluid;
  int j;
  while (grid != NULL)
  {
    if (grid->level == level && grid->cpu == CPU_Rank)
    {
      BuildCommSource_gpu(grid->fluid, fieldtype, nvar, 1);
      UPLIL(grid, nvar);
    }
    grid = grid->next;
  } /*Now all comms with source "grid" have been filled, next step is to communicate them to the pertinent destination grid*/
  FARGO_SAFE(ExecComm_gpu(ComListGhost, level, level + 1, nvar, fieldtype));
  grid = Grid_CPU_list;
  while (grid != NULL)
  {
    if (grid->level == level + 1 && grid->cpu == CPU_Rank)
    {
      BuildCommDest_gpu(grid->fluid, fieldtype, nvar, 1);
      EXECUP(grid, nvar);
    }
    grid = grid->next;
  }
  check_errors("ExecCommUpVar_gpu");
#endif
}

void ExecCommDownMeanVar_gpu(int level, int nvar, int *fieldtype)
{
#ifdef GPU
  tGrid_CPU *grid;
  grid = Grid_CPU_list;
  while (grid != NULL)
  {
    if (grid->level == level && grid->cpu == CPU_Rank)
    {
      BuildCommSource_gpu(grid->fluid, fieldtype, nvar, 1);
      DOWNMEAN(grid, nvar);
    }
    grid = grid->next;
  } /*Now all comms with source "grid" have been filled, next step is to communicate them to the pertinent destination grid*/
  FARGO_SAFE(ExecComm_gpu(ComListMean, level, level - 1, nvar, fieldtype));
  grid = Grid_CPU_list;
  while (grid != NULL)
  {
    if (grid->level == level - 1 && grid->cpu == CPU_Rank)
    {
      BuildCommDest_gpu(grid->fluid, fieldtype, nvar, 1);
      EXECDOWN(grid, nvar);
    }
    grid = grid->next;
  }
  check_errors("ExecCommDownMeanVar_gpu");
#endif
}

void ExecCommDownFlux_gpu(int level)
{
#ifdef GPU
  int nvar, fieldtype[20];
  nvar = 1 + 2 * NDIM; /* Density (or mass flux) + velocities (or momentum flux) */
#ifdef ADIABATIC
  nvar += 1; /* Energy flux */
#endif
  tGrid_CPU *grid;
  grid = Grid_CPU_list;
  FluidPatch *fluid;
  int j;
  while (grid != NULL)
  {
    if (grid->level == level && grid->cpu == CPU_Rank)
      DOWNFLUX(grid, nvar);
    grid = grid->next;
  }
  FARGO_SAFE(ExecComm_gpu(ComListFlux, level, level - 1, nvar, fieldtype));
#endif
}

void ExecComm_gpu(jCommunicator *comm, int levsrc, int levdest, int nvar, int *fieldtype)
{
#ifdef GPU
  //MPI_Request reqs[8], reqr[8];
  jCommunicator *commsave;
  MPI_Status stat;
  int n;
  int nbreq = 0;
  int parity, dimc;

  while (comm != NULL)
  {
    if ((comm->dest_level == levdest) && (comm->src_level == levsrc))
    {
      if ((comm->CPU_src != comm->CPU_dest && comm->CPU_src == CPU_Rank))
      {
#ifdef FLOAT
        MPI_Isend(comm->bufferGPU.ptr, comm->yzsize * nvar * comm->bufferGPU.pitch / sizeof(real),
                  MPI_FLOAT, comm->CPU_dest, comm->facedim,
                  DomainComm, Req + nbreq++);
#else
        MPI_Isend(comm->bufferGPU.ptr, comm->yzsize * nvar * comm->bufferGPU.pitch / sizeof(real),
                  MPI_DOUBLE, comm->CPU_dest, comm->facedim,
                  DomainComm, Req + nbreq++);
#endif
      }
      check_errors("mpi src sends to buffer");

      if ((comm->CPU_dest != comm->CPU_src && comm->CPU_dest == CPU_Rank))
      {
#ifdef FLOAT
        MPI_Irecv(comm->bufferGPU.ptr, comm->yzsize * nvar * comm->bufferGPU.pitch / sizeof(real),
                  MPI_FLOAT, comm->CPU_src, comm->facedim,
                  DomainComm, Req + nbreq++);
#else
        MPI_Irecv(comm->bufferGPU.ptr, comm->yzsize * nvar * comm->bufferGPU.pitch / sizeof(real),
                  MPI_DOUBLE, comm->CPU_src, comm->facedim,
                  DomainComm, Req + nbreq++);
#endif
      }
      check_errors("mpi dst recives buffer");
    }
    comm = comm->next;
  }
  for (n = 0; n < nbreq; n++)
    MPI_Wait(Req + n, &stat);

  MPI_Barrier(DomainComm);
  check_errors("ExecComm_gpu");
#endif
}

void ExecCommS_gpu(jCommunicator *comm, int levsrc, int levdest, int nvar, int *fieldtype)
{
#ifdef GPU
  MPI_Request reqs[8], reqr[8];
  long i, j, k, le, m, n, d;
  int field, nbreqs = 0, nbreqr = 0, staggered, stagdim, displacement[3], cut[3];
  int xmindest, ymindest, zmindest, xmaxdest, ymaxdest, zmaxdest;
  FluidPatch *fluid;
  //struct cudaPitchedPtr destiny[80];

  while (comm != NULL)
  {
    if ((comm->dest_level == levdest) && (comm->src_level == levsrc))
    {

      if (comm->CPU_src != comm->CPU_dest && comm->CPU_src == CPU_Rank)
      {
#ifdef FLOAT
        MPI_Isend(comm->bufferGPU.ptr, comm->yzsize * nvar * comm->bufferGPU.pitch / sizeof(real),
                  MPI_FLOAT, comm->CPU_dest, comm->facedim,
                  DomainComm, reqs + nbreqs++);
#else
        MPI_Isend(comm->bufferGPU.ptr, comm->yzsize * nvar * comm->bufferGPU.pitch / sizeof(real),
                  MPI_DOUBLE, comm->CPU_dest, comm->facedim,
                  DomainComm, reqs + nbreqs++);
#endif
      }
      check_errors("mpi src sends to buffer");

      if (comm->CPU_dest == CPU_Rank)
      {
        le = 0;
        xmindest = comm->imin_dest[0];
        ymindest = comm->imin_dest[1];
        zmindest = comm->imin_dest[2];
        xmaxdest = comm->imax_dest[0];
        ymaxdest = comm->imax_dest[1];
        zmaxdest = comm->imax_dest[2];
        comm->OnDst.srcPtr = comm->bufferGPU;
        fluid = comm->destg->fluid;
        BuildCommDest_gpu(fluid, fieldtype, nvar, 0);
        if (comm->CPU_dest != comm->CPU_src)
        {
#ifdef FLOAT
          MPI_Irecv(comm->bufferGPU.ptr, comm->yzsize * nvar * comm->bufferGPU.pitch / sizeof(real),
                    MPI_FLOAT, comm->CPU_src, comm->facedim,
                    DomainComm, reqr + nbreqr);
#else
          MPI_Irecv(comm->bufferGPU.ptr, comm->yzsize * nvar * comm->bufferGPU.pitch / sizeof(real),
                    MPI_DOUBLE, comm->CPU_src, comm->facedim,
                    DomainComm, reqr + nbreqr);
#endif
          MPI_Wait(reqr + nbreqr++, MPI_STATUS_IGNORE);
          /* This WAIT instruction must be here, as we have to
		 wait for the data to be ready to send it to the
		 device. We therefore split send and receive requests. */
        }
        for (le = 0; le < nvar; le++)
        {
          comm->OnDst.srcPos = make_cudaPos(0, 0, le * comm->dz);
          comm->OnDst.dstPtr = dest_gpu[le];
          cudaMemcpy3D(&(comm->OnDst));
          check_errors("cudaMemcpy3D");
        }
      }
    }
    comm = comm->next;
  }
  for (n = 0; n < nbreqs; n++)
    MPI_Wait(reqs + n, MPI_STATUS_IGNORE);
  MPI_Barrier(DomainComm);
#endif
}

void BuildCommSource_gpu(FluidPatch *fluid, int *fieldtype, int nvar, int cpy)
{
  int d, i, j, field, le = 0;
#ifdef GPU
  int centred[240];
  tGrid_CPU *grid;
  grid = fluid->desc;
  while (fluid != NULL)
  {
    if (fluid->FluidRank == Current_Fluid)
    {
      for (i = 0; i < nvar; i++)
      {
        field = fieldtype[i];
        centred[le * 3] = centred[le * 3 + 1] = centred[le * 3 + 2] = 1;
        if (field == _Density_)
        {
          JSInput_GPU(fluid->Density, __LINE__, __FILE__);
        }
        else if (field == _Energy_)
        {
          JSInput_GPU(fluid->Energy, __LINE__, __FILE__);
        }
#ifdef X
        else if (field == _Vx_)
        {
          JVInput_GPU(fluid->Velocity, _X_, __LINE__, __FILE__);
          centred[le * 3 + _X_] = 0;
        }
#endif
#ifdef Y
        else if (field == _Vy_)
        {
          JVInput_GPU(fluid->Velocity, _Y_, __LINE__, __FILE__);
          centred[le * 3 + _Y_] = 0;
        }
#endif
#ifdef Z
        else if (field == _Vz_)
        {
          JVInput_GPU(fluid->Velocity, _Z_, __LINE__, __FILE__);
          centred[le * 3 + _Z_] = 0;
        }
#endif
        source_gpu[le] = fluid->Ptr_gpu[field];   //Used for SAME
        source[le++] = fluid->Ptr_gpu[field].ptr; //Used for UP and MEAN
      }
    }
    fluid = fluid->next;
  }
  if (cpy)
  {
    cudaMemcpy(grid->centered_gpu, &centred, le * 3 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(grid->source, &source, le * sizeof(real), cudaMemcpyHostToDevice);
  }
#endif
}

void BuildCommDest_gpu(FluidPatch *fluid, int *fieldtype, int nvar, int cpy)
{
  int d, i, j, field, le = 0;
#ifdef GPU
  int centred[240];
  tGrid_CPU *grid;
  grid = fluid->desc;
  while (fluid != NULL)
  {
    if (fluid->FluidRank == Current_Fluid)
    {
      for (i = 0; i < nvar; i++)
      {
        field = fieldtype[i];
        centred[le * 3] = centred[le * 3 + 1] = centred[le * 3 + 2] = 1;
        if (field == _Density_)
        {
          JSInput_GPU(fluid->Density, __LINE__, __FILE__);
          JSOutput_GPU(fluid->Density, __LINE__, __FILE__);
        }
        else if (field == _Energy_)
        {
          JSInput_GPU(fluid->Energy, __LINE__, __FILE__);
          JSOutput_GPU(fluid->Energy, __LINE__, __FILE__);
        }
#ifdef X
        else if (field == _Vx_)
        {
          JVInput_GPU(fluid->Velocity, _X_, __LINE__, __FILE__);
          JVOutput_GPU(fluid->Velocity, _X_, __LINE__, __FILE__);
          centred[le * 3 + _X_] = 0;
        }
#endif
#ifdef Y
        else if (field == _Vy_)
        {
          JVInput_GPU(fluid->Velocity, _Y_, __LINE__, __FILE__);
          JVOutput_GPU(fluid->Velocity, _Y_, __LINE__, __FILE__);
          centred[le * 3 + _Y_] = 0;
        }
#endif
#ifdef Z
        else if (field == _Vz_)
        {
          JVInput_GPU(fluid->Velocity, _Z_, __LINE__, __FILE__);
          JVOutput_GPU(fluid->Velocity, _Z_, __LINE__, __FILE__);
          centred[le * 3 + _Z_] = 0;
        }
#endif
        dest_gpu[le] = fluid->Ptr_gpu[field];   //Used for SAME
        dest[le++] = fluid->Ptr_gpu[field].ptr; //Used for UP and MEAN
      }
    }
    fluid = fluid->next;
  }
  if (cpy)
  {
    cudaMemcpy(grid->centered_gpu, &centred, le * 3 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(grid->dest, &dest, le * sizeof(real), cudaMemcpyHostToDevice);
  }
#endif
}
