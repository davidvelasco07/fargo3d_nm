#include "fargo3d.h"


void DustDiffusion_Comm(){

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


#ifdef ALPHAVISCOSITY  
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
#endif
  
}


void DustDiffusion_Main(real dt) {

#ifndef DUSTSIZE  
  //Wait for Gas Density (QR)
  MPI_Wait(&RequestGasDensity, MPI_STATUS_IGNORE);
#ifdef ALPHAVISCOSITY
  //Wait for Gas Energy
  MPI_Wait(&RequestGasEnergy, MPI_STATUS_IGNORE);
  // This copy is needed for a local set of fluids
  MULTIFLUID( if(Fluidtype == GAS)   copy_field(Pressure, Energy);
	      if(Fluidtype == DUST)  copy_field(Energy, Pressure););
#endif
#endif

  // In principle, Diffusion_Coefficients() does not need to be called every time step
  // for temporary constant viscosity.
  FARGO_SAFE(DustDiffusion_Coefficients());

  MULTIFLUID(
	     if(Fluidtype == DUST) {
	       FARGO_SAFE(DustDiffusion_Core(dt)); // Updated dust-density is stored in Qs field.
	       FARGO_SAFE(copy_field(Density,Qs)); // Local copy of Qs to dust-density.
	       FARGO_SAFE(FillGhosts(DENS));       // Fill ghost with new density.
	       FARGO_SAFE(Reset_field(Energy));       // Fill ghost with new density.
	     });

}
