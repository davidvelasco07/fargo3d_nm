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
     

void Sources(real dt) {
     
  FARGO_SAFE(FillGhosts(PrimitiveVariables()));
  FARGO_SAFE(copy_velocities(V2VTEMP));


}

void Transport(real dt) {

  //Note: This setup does not solve any transport step.
  
  copy_velocities(VTEMP2V);
  
  return;

}
