//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Save_Face_Flux_X() {
#ifndef GPU
  INPUT(Flux);
//<EXTERNAL>
  real* flux = Flux->field_cpu;
  real* facefluxes_inf = Fluxes[_X_][0];
  real* facefluxes_sup = Fluxes[_X_][1];
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int n;
  int llINF;
  int llSUP;
//<\INTERNAL>

//<CONSTANT>
//<\CONSTANT>

//<MAIN_LOOP>

  j=k=0;
#ifdef Z  
  for(k=NGHZ; k < NGHZ+nz; k++){
#endif
#ifdef Y
    for(j=NGHY; j < NGHY+ny; j++){
#endif
//<#>
      i = NGHX + FLUXWIDTH; //the first active staggered dof is NGHX + 2
      llINF = l;                //and we want a boundary that is shared also 
                                //with the mesh at level l-1 so we go to NGHX + 3
      facefluxes_inf[(j-NGHY)+(k-NGHZ)*(ny)+FluxIndex*ny*nz] += flux[llINF]; 
      
      i = NGHX + nx - FLUXWIDTH;
      llSUP = l;
      
      facefluxes_sup[(j-NGHY)+(k-NGHZ)*(ny)+FluxIndex*ny*nz] += flux[llSUP]; 
      
//<\#>
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
  
//<\MAIN_LOOP>
#else
  Input_GPU(Flux, __LINE__, __FILE__);
  Save_Face_Flux_X_gpu();
#endif
}

void Save_Face_Flux_Y() {
#ifndef GPU
  INPUT(Flux);
//<EXTERNAL>
  real* flux = Flux->field_cpu;
  real* facefluxes_inf = Fluxes[_Y_][0];
  real* facefluxes_sup = Fluxes[_Y_][1];
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int n;
  int llINF;
  int llSUP;
//<\INTERNAL>

//<CONSTANT>
//<\CONSTANT>

//<MAIN_LOOP>

  i=k=0;
#ifdef Z  
  for(k=NGHZ; k < NGHZ+nz; k++){
#endif
#ifdef X
    for(i=NGHX; i < NGHX+nx; i++){
#endif
//<#>
      j = NGHY + FLUXWIDTH; //the first active staggered dof is NGHX + 2
      llINF = l;                //and we want a boundary that is shared also 
                                //with the mesh at level l-1 so we go to NGHX + 3
      facefluxes_inf[(i-NGHX)+(k-NGHZ)*(nx)+FluxIndex*nx*nz] += flux[llINF]; 
      
      j = NGHY + ny - FLUXWIDTH;
      llSUP = l;
      
      facefluxes_sup[(i-NGHX)+(k-NGHZ)*(nx)+FluxIndex*nx*nz] += flux[llSUP]; 
      
//<\#>
#ifdef X
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
#else
  Input_GPU(Flux, __LINE__, __FILE__);
  Save_Face_Flux_Y_gpu();
#endif
}

void Save_Face_Flux_Z() {
#ifndef GPU
  INPUT(Flux);

//<EXTERNAL>
  real* flux = Flux->field_cpu;
  real* facefluxes_inf = Fluxes[_Z_][0];
  real* facefluxes_sup = Fluxes[_Z_][1];
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
//<\EXTERNAL>

//<INTERNAL>
  int i; //Variables reserved
  int j; //for the topology
  int k; //of the kernels
  int n;
  int llINF;
  int llSUP;
//<\INTERNAL>

//<CONSTANT>
//<\CONSTANT>

//<MAIN_LOOP>

  i=j=0;
#ifdef Y  
  for(j=NGHY; j < NGHY+ny; j++){
#endif
#ifdef X
    for(i=NGHX; i < NGHX+nx; i++){
#endif
//<#>
      k = NGHZ + FLUXWIDTH;
      llINF = l;           

      facefluxes_inf[(i-NGHX)+(j-NGHY)*(nx)+FluxIndex*nx*ny] += flux[llINF]; 
      
      k = NGHZ + nz - FLUXWIDTH;
      llSUP = l;
      
      facefluxes_sup[(i-NGHX)+(j-NGHY)*(nx)+FluxIndex*nx*ny] += flux[llSUP]; 
      
//<\#>
#ifdef X
    }
#endif
#ifdef Y
  }
#endif
//<\MAIN_LOOP>  
#else
  Input_GPU(Flux, __LINE__, __FILE__);
  Save_Face_Flux_Z_gpu();
#endif
}

