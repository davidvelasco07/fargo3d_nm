//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void _collisions_cpu(real dt, int id1, int id2, int id3, int option) {
};

//Prototypes -------------------------
void ComputeCBcollisions_c_cpu(real);
void ComputeCBcollisions_cx_cpu(real);
void ComputeCBcollisions_cy_cpu(real);
void ComputeCBcollisions_cz_cpu(real);
void UpdateVelcollisions_x_cpu(real);
void UpdateVelcollisions_y_cpu(real);
void UpdateVelcollisions_z_cpu(real);
void ComputeDragCoeff_cpu();
// -----------------------------------


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
  MPI_Bcast(Total_Density->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ), MPI_DOUBLE, 0, FluidsComm);

  // Centered drag coefficient pre-factor
  FARGO_SAFE(ComputeDragCoeff_cpu());
  //--------------------------------------------------------------------------------------
  
  Reset_field_cpu (Slope);
  //First step: Compute the constant C and B for each fluids  
  MULTIFLUID(ComputeCBcollisions_c_cpu(dt)); // local C and local B (x, y, z)
  MPI_Iallreduce(MPI_IN_PLACE, Slope->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_c); // obtain C

  //Second step: Compute the total C and B values (sum for all fluids)
#ifdef X
  Reset_field_cpu (Mpx);
  MULTIFLUID(ComputeCBcollisions_cx_cpu(dt)); // local B (x)
  MPI_Iallreduce(MPI_IN_PLACE, Mpx->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_x); // obtain B (x)
#endif
#ifdef Y  
  Reset_field_cpu (Mpy);
  MULTIFLUID(ComputeCBcollisions_cy_cpu(dt)); // local B (y)
  MPI_Iallreduce(MPI_IN_PLACE, Mpy->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_y); // obtain B (y)  
#endif
#ifdef Z  
  Reset_field_cpu (Mpz);
  MULTIFLUID(ComputeCBcollisions_cz_cpu(dt)); // local B (z)
  MPI_Iallreduce(MPI_IN_PLACE, Mpz->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_z); // obtain B (z)    
#endif
  
  MPI_Wait(&request_c, MPI_STATUS_IGNORE);
  
  //Third step: Update the velocities of each fluid
#ifdef X
  MPI_Wait(&request_x, MPI_STATUS_IGNORE);
  MULTIFLUID(UpdateVelcollisions_x_cpu(dt)); //update velocities
#endif

#ifdef Y
  MPI_Wait(&request_y, MPI_STATUS_IGNORE);
  MULTIFLUID(UpdateVelcollisions_y_cpu(dt)); //update velocities
#endif
  
#ifdef Z  
  MPI_Wait(&request_z, MPI_STATUS_IGNORE);
  MULTIFLUID(UpdateVelcollisions_z_cpu(dt)); //update velocities
#endif
  
}
