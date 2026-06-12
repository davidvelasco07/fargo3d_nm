#include "fargo3d.h"

MPI_Comm RadialCPUBeam;
boolean  CommCreated = NO;

void CreateRadialComm() {
  if (CommCreated == NO) {
    MPI_Comm_split (MPI_COMM_WORLD, K, J, &RadialCPUBeam);
    CommCreated = YES;
  }
}

void InterCPUScan_cpu () {
  CreateRadialComm ();
  Input2D_CPU (LocalDepth, __LINE__, __FILE__);
  Output2D_CPU (LocalDepth, __LINE__, __FILE__);
#ifdef FLOAT
  MPI_Scan (MPI_IN_PLACE, LocalDepth->field_cpu, (Nx+2*NGHX)*(Nz+2*NGHZ), MPI_FLOAT, MPI_SUM, RadialCPUBeam);
#else
  MPI_Scan (MPI_IN_PLACE, LocalDepth->field_cpu, (Nx+2*NGHX)*(Nz+2*NGHZ), MPI_DOUBLE, MPI_SUM, RadialCPUBeam);
#endif
}

void InterCPUScan_gpu () {
  /* This function should not be invoked yet, as present
     implementations of CUDA aware MPI do not support reduction type
     operations (MPI_Reduce, MPI_Allreduce, MPI_Scan). 'change_arch.c'
     has been edited to take this limitation into account. */

  CreateRadialComm ();
  Input2D_GPU (LocalDepth, __LINE__, __FILE__);
  Output2D_GPU (LocalDepth, __LINE__, __FILE__);
#ifdef FLOAT
  MPI_Scan (MPI_IN_PLACE, LocalDepth->field_gpu, (LocalDepth->pitch)*(Nz+2*NGHZ), MPI_FLOAT, MPI_SUM, RadialCPUBeam);
#else
  MPI_Scan (MPI_IN_PLACE, LocalDepth->field_gpu, (LocalDepth->pitch)*(Nz+2*NGHZ), MPI_DOUBLE, MPI_SUM, RadialCPUBeam);
#endif
}
