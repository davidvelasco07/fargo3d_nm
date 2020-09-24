//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void _collisions_cpu(real dt, int id1, int id2, int id3, int option) {
};

void Collisions(real dt, int option) {

  MPI_Request request_c;
  MPI_Request request_x;
  MPI_Request request_y;
  MPI_Request request_z;

  //--------------------------------------------------------------------------------------
  /* 
     Note: We should try to send the gas density at the very begining
     of the time step so we hide its communication time.
  */
 
  MULTIFLUID(if(Fluidtype == GAS) copy_field(Total_Density, Density));
#ifdef MPICUDA
  MPI_Bcast(Total_Density->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ), MPI_DOUBLE, 0, FluidsComm);
#else
  INPUT(Total_Density);
  MPI_Bcast(Total_Density->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ), MPI_DOUBLE, 0, FluidsComm);
#endif
  // Centered drag coefficient pre-factor
  FARGO_SAFE(ComputeDragCoeff());
  //--------------------------------------------------------------------------------------

  Reset_field(Slope);
  //First step: Compute the constant C and B for each fluids  
  MULTIFLUID(ComputeCBcollisions_c(dt)); // local C and local B (x, y, z)
#ifdef MPICUDA
  MPI_Iallreduce(MPI_IN_PLACE, Slope->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_c); // obtain C
#else
  INPUT(Slope);
  MPI_Iallreduce(MPI_IN_PLACE, Slope->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_c); // obtain C
#endif

  //Second step: Compute the total C and B values (sum for all fluids)
#ifdef X
  Reset_field(Mpx);
  MULTIFLUID(ComputeCBcollisions_cv(dt,0)); // local B (x)
#ifdef MPICUDA
  MPI_Iallreduce(MPI_IN_PLACE, Mpx->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_x); // obtain B (x)
#else

  INPUT(Mpx);
  MPI_Iallreduce(MPI_IN_PLACE, Mpx->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_x); // obtain B (x)
#endif
#endif
#ifdef Y  
  Reset_field(Mpy);
  MULTIFLUID(ComputeCBcollisions_cv(dt,1)); // local B (y)
#ifdef MPICUDA
  MPI_Iallreduce(MPI_IN_PLACE, Mpy->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_y); // obtain B (y)
#else
  INPUT(Mpy);
  MPI_Iallreduce(MPI_IN_PLACE, Mpy->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_y); // obtain B (y)  
#endif
#endif
#ifdef Z  
  Reset_field_cpu (Mpz);
  MULTIFLUID(ComputeCBcollisions_cv(dt,2)); // local B (z)
#ifdef MPICUDA
  MPI_Iallreduce(MPI_IN_PLACE, Mpz->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_z); // obtain B (z)
#else
  INPUT(Mpz);
  MPI_Iallreduce(MPI_IN_PLACE, Mpz->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_z); // obtain B (z)    
#endif
#endif
 
  MPI_Wait(&request_c, MPI_STATUS_IGNORE);
  //Third step: Update the velocities of each fluid
#ifdef X
  MPI_Wait(&request_x, MPI_STATUS_IGNORE);
  MULTIFLUID(UpdateVelcollisions(dt,0)); //update velocities
#endif

#ifdef Y
  MPI_Wait(&request_y, MPI_STATUS_IGNORE);
  MULTIFLUID(UpdateVelcollisions(dt,1)); //update velocities
#endif
  
#ifdef Z  
  MPI_Wait(&request_z, MPI_STATUS_IGNORE);
  MULTIFLUID(UpdateVelcollisions(dt,2)); //update velocities
#endif

}
