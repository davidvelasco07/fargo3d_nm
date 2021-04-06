#include "fargo3d.h"


void DustDiffusion_Comm(){

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


void DustDiffusion_Main(real dt) {

#ifndef DUSTSIZE  
  //Wait for Gas Density (QR)
  MPI_Wait(&RequestGasDensity, MPI_STATUS_IGNORE);
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
