#include "fargo3d.h"

static void (*__VanLeerX) (Field *, Field *, Field *, real);

void VanLeerX(Field *Q, Field *Qs, Field *Vx_t, real dt) {
  FARGO_SAFE(VanLeerX_a(Q));
  FARGO_SAFE(VanLeerX_b(dt, Q, Qs, Vx_t));
}


void TransportX(Field *Q, Field *Qs, Field *Vx_t, real dt) { 
  FARGO_SAFE(DivideByRho(Q,DivRho));
  __VanLeerX(DivRho, Qs, Vx_t, dt);
  FARGO_SAFE(Update_a_X (dt,Qs, Vx_t, Flux));
  FARGO_SAFE(Save_Face_Flux_X ()); // <== will be used by coarser level
  FARGO_SAFE(OverWriteBoundaryFluxes (_X_)); // <== we use information from finer level
  FARGO_SAFE(Update_b_X (Q, Flux));
}
void TransportY(Field *Q, Field *Qs, real dt) {
  FARGO_SAFE(DivideByRho(Q,DivRho));
  FARGO_SAFE(VanLeerY_a(DivRho));
  FARGO_SAFE(VanLeerY_b(dt, DivRho, Qs));
  FARGO_SAFE(Update_a_Y (dt, Qs, Flux));
  FARGO_SAFE(Save_Face_Flux_Y ()); // <== will be used by coarser level
  FARGO_SAFE(OverWriteBoundaryFluxes (_Y_)); // <== we use information from finer level
  FARGO_SAFE(Update_b_Y (Q, Flux));
}

void TransportZ(Field *Q, Field *Qs, real dt) {
  FARGO_SAFE(DivideByRho(Q,DivRho));
  FARGO_SAFE(VanLeerZ_a(DivRho));
  FARGO_SAFE(VanLeerZ_b(dt, DivRho, Qs));
  FARGO_SAFE(Update_a_Z (dt, Qs, Flux));
  FARGO_SAFE(Save_Face_Flux_Z ()); // <== will be used by coarser level
  FARGO_SAFE(OverWriteBoundaryFluxes (_Z_)); // <== we use information from finer level
  FARGO_SAFE(Update_b_Z (Q, Flux));
}



void X_advection (Field *Vx_t, real dt) {
#ifdef X
  FluxIndex = 0;
  __VanLeerX(Density, DensStar, Vx_t, dt);
  TransportX(Mpx, Qs, Vx_t, dt);
  TransportX(Mmx, Qs, Vx_t, dt);
#endif
#ifdef Y
  TransportX(Mpy, Qs, Vx_t, dt);
  TransportX(Mmy, Qs, Vx_t, dt);
#endif
#ifdef Z
  TransportX(Mpz, Qs, Vx_t, dt);
  TransportX(Mmz, Qs, Vx_t, dt);
#endif
#ifdef ADIABATIC
  TransportX(Energy, Qs, Vx_t, dt);
#endif
#ifdef LABELED
  TransportX(Label, Qs, Vx_t, dt);
#endif
  Update_dens_a_X (dt, DensStar, Vx_t, Flux);
  Save_Face_Flux_X (); 
  OverWriteBoundaryFluxes (_X_);
  Update_b_X (Density, Flux);
}

void transport(real dt){
 
#ifdef X
  FARGO_SAFE(momenta_x());
#endif
#ifdef Y
  FARGO_SAFE(momenta_y());
#endif
#ifdef Z
  FARGO_SAFE(momenta_z());
#endif
  
#ifdef Z
  FARGO_SAFE(VanLeerZ_a(Density));
  FARGO_SAFE(VanLeerZ_b(dt, Density, DensStar));
  FluxIndex = 0;
#ifdef X
  TransportZ(Mpx, Qs, dt);
  TransportZ(Mmx, Qs, dt);
#endif
#ifdef Y
  TransportZ(Mpy, Qs, dt);
  TransportZ(Mmy, Qs, dt);
#endif
#ifdef Z
  TransportZ(Mpz, Qs, dt);
  TransportZ(Mmz, Qs, dt);
#endif
#ifdef ADIABATIC
  TransportZ(Energy, Qs, dt);
#endif
#ifdef LABELED
  TransportZ(Label, Qs, dt);
#endif
  Update_dens_a_Z (dt, DensStar, Flux);
  Save_Face_Flux_Z (); // <== will be used by coarser level
  OverWriteBoundaryFluxes (_Z_); // <== we use information from finer level
  Update_b_Z (Density, Flux);
#endif
 
#ifdef Y
  FARGO_SAFE(VanLeerY_a(Density));
  FARGO_SAFE(VanLeerY_b(dt, Density, DensStar));
  FluxIndex = 0;
#ifdef X  
  TransportY(Mpx, Qs, dt);
  TransportY(Mmx, Qs, dt);
#endif
#ifdef Y
  TransportY(Mpy, Qs, dt);
  TransportY(Mmy, Qs, dt);
#endif
#ifdef Z
  TransportY(Mpz, Qs, dt);
  TransportY(Mmz, Qs, dt);
#endif
#ifdef ADIABATIC
  TransportY(Energy, Qs, dt);
#endif
#ifdef LABELED
  TransportY(Label, Qs, dt);
#endif
  Update_dens_a_Y (dt, DensStar, Flux);
  Save_Face_Flux_Y ();
  OverWriteBoundaryFluxes (_Y_); 
  Update_b_Y (Density, Flux);
#endif
 
#ifdef X
#ifdef STANDARD
  __VanLeerX = VanLeerX;
  X_advection (Vx_temp, dt);
#else
  if (Current_Level > 0) {
    __VanLeerX = VanLeerX;
    X_advection (Vx_temp, dt);
  } else {
    FARGO_SAFE(ComputeResidual(dt));
    __VanLeerX = VanLeerX;
    FirstPassXAdvection = YES;
    X_advection (Vx, dt); // Vx => variable residual
    FirstPassXAdvection = NO;
    X_advection (Vx_temp, dt); // Vx_temp => fixed residual @ given r.
    AdvectSHIFT(Mpx, Nshift);
    AdvectSHIFT(Mmx, Nshift);
#ifdef Y
    AdvectSHIFT(Mpy, Nshift);
    AdvectSHIFT(Mmy, Nshift);
#endif
#ifdef Z
    AdvectSHIFT(Mpz, Nshift);
    AdvectSHIFT(Mmz, Nshift);
#endif
#ifdef ADIABATIC
    AdvectSHIFT(Energy, Nshift);
#endif
#ifdef LABELED
    AdvectSHIFT(Label, Nshift);
#endif
    AdvectSHIFT(Density, Nshift);
  }
#endif
#endif
  FARGO_SAFE(NewVelocity_x());
  FARGO_SAFE(NewVelocity_y());
  FARGO_SAFE(NewVelocity_z());
}
