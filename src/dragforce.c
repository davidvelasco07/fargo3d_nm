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
  MPI_Iallreduce(QL->field_gpu, QR->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasDensity); 
#else
  INPUT(QL);
  OUTPUT(QR);
  MPI_Iallreduce(QL->field_cpu, QR->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
  		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasDensity); 
#endif
  
  //Send Gas Energy to Dust-Fluids
#ifdef MPICUDA
  MPI_Iallreduce(MPI_IN_PLACE, Energy->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasEnergy); 
#else
  INPUT(Energy);
  OUTPUT(Energy);
  MPI_Iallreduce(MPI_IN_PLACE, Energy->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
  		 MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasEnergy);
#endif

  
}

void DragForce_a(real dt) {
  
#ifdef DUSTSIZE  
  //Wait for Gas Density (QR)
  MPI_Wait(&RequestGasDensity, MPI_STATUS_IGNORE);
  //Wait for Gas Energy  (DivRho)
  MPI_Wait(&RequestGasEnergy, MPI_STATUS_IGNORE);
  // This copy is needed for a local set of fluids
  MULTIFLUID( if(Fluidtype == GAS)   copy_field(DivRho, Energy);
	      if(Fluidtype == DUST)  copy_field(Energy, DivRho););
#endif
  
  // Centered drag coefficient pre-factor
  FARGO_SAFE(DragForce_Coeff()); //store dragcoeff in Qs
  //-------------------------------------------------------------------------------------

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
}

void DragForce_b(real dt) {

  MPI_Request request_x;
  MPI_Request request_y;
  MPI_Request request_z;
  
  //Second step: 
  // Compute local coefficient CV_i (Mmx,Mmy,Mmz) and obtain global sum CV (Mpx,Mpy,Mpz)
#ifdef X
  Reset_field(Mmx);
  MULTIFLUID(DragForce_SumCV(dt,0)); 
#ifdef MPICUDA
  MPI_Iallreduce(Mmx->field_gpu, Mpx->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_x); 
#else

  INPUT(Mmx);
  OUTPUT(Mpx);
  MPI_Iallreduce(Mmx->field_cpu, Mpx->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_x); 
#endif
#endif

#ifdef Y  
  Reset_field(Mmy);
  MULTIFLUID(DragForce_SumCV(dt,1)); 
#ifdef MPICUDA
  MPI_Iallreduce(Mmy->field_gpu, Mpy->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_y); 
#else
  INPUT(Mmy);
  OUTPUT(Mpy);
  MPI_Iallreduce(Mmy->field_cpu, Mpy->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_y);  
#endif
#endif
#ifdef Z  
  Reset_field(Mmz); 
  MULTIFLUID(DragForce_SumCV(dt,2)); 
#ifdef MPICUDA
  MPI_Iallreduce(Mmz->field_gpu, Mpz->field_gpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		 MPI_DOUBLE, MPI_SUM, FluidsComm, &request_z); 
#else
  INPUT(Mpz);
  OUTPUT(Mmz);
  MPI_Iallreduce(Mmz->field_cpu, Mpz->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
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
