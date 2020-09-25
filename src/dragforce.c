//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void DragForce(real dt) {

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
  FARGO_SAFE(DragForce_Coeff());
  //--------------------------------------------------------------------------------------

  Reset_field(Slope);  
  MULTIFLUID(DragForce_SumC(dt)); 
#ifdef MPICUDA
  MPI_Iallreduce(MPI_IN_PLACE, Slope->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_c); // obtain C
#else
  INPUT(Slope);
  MPI_Iallreduce(MPI_IN_PLACE, Slope->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_c); // obtain C
#endif

  //Second step: Compute the total CV (sum for all fluids)
#ifdef X
  Reset_field(Mpx);
  MULTIFLUID(DragForce_SumCV(dt,0)); 
#ifdef MPICUDA
  MPI_Iallreduce(MPI_IN_PLACE, Mpx->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_x); 
#else

  INPUT(Mpx);
  MPI_Iallreduce(MPI_IN_PLACE, Mpx->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_x); 
#endif
#endif
#ifdef Y  
  Reset_field(Mpy);
  MULTIFLUID(DragForce_SumCV(dt,1)); 
#ifdef MPICUDA
  MPI_Iallreduce(MPI_IN_PLACE, Mpy->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_y); 
#else
  INPUT(Mpy);
  MPI_Iallreduce(MPI_IN_PLACE, Mpy->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_y);  
#endif
#endif
#ifdef Z  
  Reset_field_cpu (Mpz);
  MULTIFLUID(DragForce_SumCV(dt,2)); 
#ifdef MPICUDA
  MPI_Iallreduce(MPI_IN_PLACE, Mpz->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_z); 
#else
  INPUT(Mpz);
  MPI_Iallreduce(MPI_IN_PLACE, Mpz->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_z);     
#endif
#endif
 
  MPI_Wait(&request_c, MPI_STATUS_IGNORE);
  //Third step: Update the velocities of each fluid
#ifdef X
  MPI_Wait(&request_x, MPI_STATUS_IGNORE);
  MULTIFLUID(DragForce_UpdateVel(dt,0)); //update velocities
#endif

#ifdef Y
  MPI_Wait(&request_y, MPI_STATUS_IGNORE);
  MULTIFLUID(DragForce_UpdateVel(dt,1)); //update velocities
#endif
  
#ifdef Z  
  MPI_Wait(&request_z, MPI_STATUS_IGNORE);
  MULTIFLUID(DragForce_UpdateVel(dt,2)); //update velocities
#endif

}
