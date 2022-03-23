#include "fargo3d.h"


void DustDiffusion_Comm(){

  if(NFLUIDCOLORS>1){
    Reset_field(QL);
    #ifdef ADIABATIC
    Reset_field(QLE);
    #endif
    MULTIFLUID( if(Fluidtype == GAS){
      copy_field(QL , Density);
      #ifdef ADIABATIC
      copy_field(QLE, Energy);
      #endif
      });
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
    MPI_Iallreduce(QL->field_cpu, QR->field_cpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
    MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasDensity); 
    #ifdef ADIABATIC
    INPUT(QLE);
    OUTPUT(QRE);
    MPI_Iallreduce(QLE->field_cpu, QRE->field_cpu, (Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ),
    MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestGasEnergy); 
    #endif
    #endif
  }
  else{
    MULTIFLUID( if(Fluidtype == GAS){
    copy_field(QR , Density);
    #ifdef ADIABATIC
    copy_field(QRE, Energy);
    #endif
    });
  }
}


void DustDiffusion_Main(real dt) {
#ifndef DUSTSIZE
  if(NFLUIDCOLORS>1){  
  //Wait for Gas Density (QR)
  MPI_Wait(&RequestGasDensity, MPI_STATUS_IGNORE);
  #ifdef ADIABATIC
  //Wait for Gas Energy (QRE)
  MPI_Wait(&RequestGasEnergy, MPI_STATUS_IGNORE);
  #endif
  }
#endif
  // In principle, Diffusion_Coefficients() does not need to be called every time step
  // for temporary constant viscosity.
  FARGO_SAFE(DustDiffusion_Coefficients());

  MULTIFLUID(
	     if(Fluidtype == DUST) {
	       FARGO_SAFE(DustDiffusion_Core(dt)); // Updated dust-density is stored in Qs field.
	       FARGO_SAFE(copy_field(Density,Qs)); // Local copy of Qs to dust-density.
	       FARGO_SAFE(FillGhosts(DENS));       // Fill ghost with new density.
	     });
}
