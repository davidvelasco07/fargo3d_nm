#include "fargo3d.h"

TimeProcess t_Comm;
TimeProcess t_Hydro;
TimeProcess t_Mhd;
TimeProcess t_sub1;
TimeProcess t_sub1_x;
TimeProcess t_sub1_y;
TimeProcess t_sub1_z;

void FillGhosts (int var) {

  InitSpecificTime (&t_Comm, "MPI Communications");
  FARGO_SAFE(comm (var));
  GiveSpecificTime (t_Comm);
  FARGO_SAFE(boundaries()); // Always after a comm.

#if defined(Y)
  if (NY == 1)    /* Y dimension is mute */
    CheckMuteY();
#endif
#if defined(Z)
  if (NZ == 1)    /* Z dimension is mute */
    CheckMuteZ();
#endif

}

static boolean Resistivity_Profiles_Filled = NO;

void Fill_Resistivity_Profiles () {

  OUTPUT2D(Eta_profile_xi);
  OUTPUT2D(Eta_profile_xizi);
  OUTPUT2D(Eta_profile_zi);

  int j,k;
  if (Resistivity_Profiles_Filled) return;
  real* eta_profile_xi = Eta_profile_xi->field_cpu;
  real* eta_profile_xizi = Eta_profile_xizi->field_cpu;
  real* eta_profile_zi = Eta_profile_zi->field_cpu;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      eta_profile_xi[l2D] = Resistivity (Ymin(j),Zmed(k));
      eta_profile_xizi[l2D] = Resistivity (Ymin(j),Zmin(k));
      eta_profile_zi[l2D] = Resistivity (Ymed(j),Zmin(k));
    }
  }
  Resistivity_Profiles_Filled = YES;
}


void AlgoGas() {
  
  real dtemp=0.0;
  real dt=1.0;  
  int var=0;

  while(dtemp<DT) { // DT LOOP

    dt = (1.0)/(real)NX*.5;

    FARGO_SAFE(copy_velocities(V2VTEMP));
    transport(dt);

    dtemp += dt;

    if(CPU_Master) {
      if (FullArrayComms)
	printf("%s", "!");
      else {
	if (ContourComms)
	  printf("%s", ":");
	else
	  printf("%s", ".");
      }
#ifndef NOFLUSH
      fflush(stdout);
#endif
    }
    if (ForwardOneStep == YES) prs_exit(EXIT_SUCCESS);
    PhysicalTime+=dt;
    FullArrayComms = 0;
    ContourComms = 0;

    FARGO_SAFE(FillGhosts (PrimitiveVariables()));


  }

  dtemp = 0.0;
  if(CPU_Master) printf("%s", "\n");

}
