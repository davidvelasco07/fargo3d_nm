#include "fargo3d.h"

void RT_message () {
  static boolean FirstTime = YES;
  if (FirstTime) {
    FirstTime = NO;
    masterprint ("\n\n\n////////////////////////\n");
    masterprint ("Radiative Tranfer Module\n");
    masterprint ("////////////////////////\n\n");
    masterprint ("It has only been tested in spherical\n");
    masterprint ("coordinates at the present time\n");
#ifdef IRRADIATION
    masterprint ("\n\n\n////////////////////////\n");
    masterprint ("Irradiation submodule activated\n");
    masterprint ("////////////////////////\n\n");
    masterprint ("It has also been tested only in spherical\n");
    masterprint ("coordinates at the present time\n");
    masterprint ("It assumes a point-like radiation source at the mesh center\n");
#endif    
  }
}

void RT_main (real dt) {
  RT_message ();
#ifdef IRRADIATION
  FARGO_SAFE (RT_IRR_TauStellar ());
  FARGO_SAFE (ComputeStellarHeating());
#endif
  FARGO_SAFE (ComputeTemperature());
  FARGO_SAFE (ComputeOpacity());
  FARGO_SAFE (Heating()); //Viscous and stellar heating into Qplus.
  ExecCommSame (Current_Level, QPLUS);
  // Here we need to sync the buffer zones of Qplus
  FARGO_SAFE (ComputeEtas (dt));
  FARGO_SAFE (ComputeDiffusionCoefficients ());
  FARGO_SAFE (ComputeMatrixElements (dt));
  FARGO_SAFE (SetRTBoundaryConditions ());/*?*/
  //FARGO_SAFE (PredictNewEnergyRad (dt));
  FARGO_SAFE (SolveMatrix ());
  //FARGO_SAFE (NewDiffEnergyRad (1.0/dt));
  FARGO_SAFE (RadEnergyToTemperature ());
  FARGO_SAFE (TemperatureToEnergy());
  /* The RT module has only updated the energy. All other primitive
     variables are left untouched. We therefore only need to
     synchronize the internal energy across CPUs */
}
