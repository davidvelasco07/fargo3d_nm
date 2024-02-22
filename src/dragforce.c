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
  #ifdef ADIABATIC
  Reset_field(QLE);
  #endif
  MULTIFLUID( if(Fluidtype == GAS){
    copy_field(QL , Density);
    #ifdef ADIABATIC
    copy_field(QLE, Energy);
    #endif
    }
    );
  
#ifdef MPICUDA
  MPI_Iallreduce(QL->field_gpu, QR->field_gpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasDensity); 
#ifdef ADIABATIC
  MPI_Iallreduce(QLE->field_gpu, QRE->field_gpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasEnergy); 
#endif
#else
  INPUT(QL);
  OUTPUT(QR);
  #ifdef ADIABATIC
  INPUT(QLE);
  OUTPUT(QRE);
  #endif
  MPI_Iallreduce(QL->field_cpu, QR->field_cpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
  		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasDensity); 
  #ifdef ADIABATIC
  MPI_Iallreduce(QLE->field_cpu, QRE->field_cpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasEnergy); 
  #endif
#endif
  
}

void DragForce_a(real dt){}

void DragForce_b(real dt) {

  MPI_Request request_c;
  MPI_Request request_x;
  MPI_Request request_y;
  MPI_Request request_z;

     
#ifdef DUSTSIZE
  //Wait for Gas Density (QR)
  MPI_Wait(&RequestGasDensity, MPI_STATUS_IGNORE);
#ifdef ADIABATIC
//Wait for Gas Energy (QRE)
  MPI_Wait(&RequestGasEnergy, MPI_STATUS_IGNORE);
#endif
#endif

// Centered drag coefficient pre-factor
FARGO_SAFE(DragForce_Coeff()); //store dragcoeff in Qs
//-------------------------------------------------------------------------------------

#ifdef FEEDBACK
  //Compute local coefficient C_i (DensStar) and obtain global sum C (slope)
  Reset_field(DensStar);  
  MULTIFLUID(DragForce_SumC(dt)); 
#ifdef MPICUDA
  MPI_Iallreduce(DensStar->field_gpu, Slope->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_c); // obtain C
#else
  INPUT(DensStar);
  OUTPUT(Slope);
  MPI_Iallreduce(DensStar->field_cpu, Slope->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_c); // obtain C
#endif
#endif

#ifdef X
  Reset_field(Mmx);
#ifdef FEEDBACK 
  MULTIFLUID(DragForce_SumCV(dt,0)); 
#else
  MULTIFLUID( if(Fluidtype == GAS) copy_field(Mmx, Vx_temp) );
#endif
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
#ifdef FEEDBACK 
  MULTIFLUID(DragForce_SumCV(dt,1)); 
#else
  MULTIFLUID( if(Fluidtype == GAS) copy_field(Mmy, Vy_temp) );
#endif
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
#ifdef FEEDBACK 
  MULTIFLUID(DragForce_SumCV(dt,2)); 
#else
  MULTIFLUID( if(Fluidtype == GAS) copy_field(Mmz, Vz_temp) );
#endif
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

#ifdef FEEDBACK
  MPI_Wait(&request_c, MPI_STATUS_IGNORE);
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
