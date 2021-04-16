//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void DragForce_Comm(){

  //Communicate gas density and energy for the drag coefficient calculation
  //Send Gas Density to Dust-Fluids
  Reset_field(QL);
  MULTIFLUID( if(Fluidtype == GAS) copy_field(QL, Density) );
  
#ifdef MPICUDA
  MPI_Iallreduce(QL->field_gpu, QR->field_gpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasDensity); 
#else
  INPUT(QL);
  OUTPUT(QR);
  MPI_Iallreduce(QL->field_cpu, QR->field_cpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
  		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasDensity); 
#endif
  
}

void DragForce_a(real dt) {
  
}

void DragForce_b(real dt) {

  MPI_Request request_x;
  MPI_Request request_y;
  MPI_Request request_z;

#ifdef DUSTSIZE
  //Wait for Gas Density (QR)
  MPI_Wait(&RequestGasDensity, MPI_STATUS_IGNORE);
#endif

  // Centered drag coefficient pre-factor
  FARGO_SAFE(DragForce_Coeff()); //store dragcoeff in Qs
  //-------------------------------------------------------------------------------------

    
#ifdef X
  Reset_field(Mmx);
  MULTIFLUID( if(Fluidtype == GAS) copy_field(Mmx, Vx_temp) );

#ifdef MPICUDA
  MPI_Iallreduce(Mmx->field_gpu, Mpx->field_gpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_x); 
#else
  INPUT(Mmx);
  OUTPUT(Mpx);
  MPI_Iallreduce(Mmx->field_cpu, Mpx->field_cpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_x); 
#endif
#endif
  
#ifdef Y  
  Reset_field(Mmy);
  MULTIFLUID( if(Fluidtype == GAS) copy_field(Mmy, Vy_temp) );
#ifdef MPICUDA
  MPI_Iallreduce(Mmy->field_gpu, Mpy->field_gpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_y); 
#else
  INPUT(Mmy);
  OUTPUT(Mpy);
  MPI_Iallreduce(Mmy->field_cpu, Mpy->field_cpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_y);  
#endif
#endif
#ifdef Z  
  Reset_field(Mmz); 
  MULTIFLUID( if(Fluidtype == GAS) copy_field(Mmz, Vz_temp) );
#ifdef MPICUDA
  MPI_Iallreduce(Mmz->field_gpu, Mpz->field_gpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_z); 
#else
  INPUT(Mmz);
  OUTPUT(Mpz);
  MPI_Iallreduce(Mmz->field_cpu, Mpz->field_cpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_z);     
#endif
#endif


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

  MULTIFLUID(copy_velocities(VTEMP2V));  
}
