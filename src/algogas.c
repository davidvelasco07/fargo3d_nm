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
  FARGO_SAFE(ExecCommSame (Current_Level, var));
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
     
  SetupHook1 (); //Setup specific hook. Defaults to empty function.
  
    
  InitSpecificTime (&t_Hydro, "Eulerian Hydro (no transport) algorithms");
  
  // REGARDLESS OF WHETHER WE USE FARGO, Vx IS ALWAYS THE TOTAL VELOCITY IN X

#ifdef POTENTIAL

  if(FluidIndex==0){
    //We compute the potential only once per timestep
    FARGO_SAFE(compute_potential(dt));
  }
  else{
    if (Corotating) {
      //We correct for the rest of the fluids
      FARGO_SAFE(CorrectVtheta(Domega));
    }
  }
#endif

  //Equations of state-----------------------------------------------------------
  if (Fluidtype == GAS){
#ifdef ADIABATIC
    FARGO_SAFE(ComputePressureFieldAd());
#endif
#ifdef ISOTHERMAL
    FARGO_SAFE(ComputePressureFieldIso());
#endif
#ifdef POLYTROPIC
    FARGO_SAFE(ComputePressureFieldPoly());
#endif
  }
  if(Fluidtype==DUST) FARGO_SAFE(Reset_field(Pressure));
  //-----------------------------------------------------------------------------

  
#if ((defined(SHEARINGSHEET2D) || defined(SHEARINGBOX3D)) && !defined(SHEARINGBC))
  FARGO_SAFE(NonReflectingBC(Vy));
#endif

#ifdef X
  FARGO_SAFE(SubStep1_x(dt));
#endif    
#ifdef Y
  FARGO_SAFE(SubStep1_y(dt));
#endif  
#ifdef Z
  FARGO_SAFE(SubStep1_z(dt));
#endif
  
#if (defined(VISCOSITY) || defined(ALPHAVISCOSITY))
  if (Fluidtype == GAS) viscosity(dt);
#endif
  
#ifndef NOSUBSTEP2
  FARGO_SAFE(SubStep2_a(dt));
  FARGO_SAFE(SubStep2_b(dt));
#endif

  // NOW: Vx INITIAL X VELOCITY, Vx_temp UPDATED X VELOCITY FROM SOURCE TERMS + ARTIFICIAL VISCOSITY

#ifdef ADIABATIC
  if (Fluidtype == GAS)
    FARGO_SAFE(SubStep3(dt));
#endif
  FARGO_SAFE(copy_velocities(VTEMP2V));  
  GiveSpecificTime (t_Hydro);
  
#ifdef MHD //-------------------------------------------------------------------
  if(Fluidtype == GAS){
    InitSpecificTime (&t_Mhd, "MHD algorithms");
    FARGO_SAFE(copy_velocities(VTEMP2V));
#ifndef STANDARD // WE USE THE FARGO ALGORITHM
    FARGO_SAFE(ComputeVmed(Vx));
    FARGO_SAFE(ChangeFrame(-1, Vx, VxMed)); //Vx becomes the residual velocity
    VxIsResidual = YES;
#endif
     
    ComputeMHD(dt);

#ifndef STANDARD
    FARGO_SAFE(ChangeFrame(+1, Vx, VxMed)); //Vx becomes the total, updated velocity
    VxIsResidual = NO;
#endif //STANDARD
    FARGO_SAFE(copy_velocities(V2VTEMP));
    // THIS COPIES Vx INTO Vx_temp
    GiveSpecificTime (t_Mhd);
  }
#endif //END MHD----------------------------------------------------------------

  InitSpecificTime (&t_Hydro, "Transport algorithms");

#if ((defined(SHEARINGSHEET2D) || defined(SHEARINGBOX3D)) && !defined(SHEARINGBC))
  FARGO_SAFE(NonReflectingBC (Vy_temp));
#endif

#ifdef MHD //-------------------------------------------------------------------
  if(Fluidtype == GAS){ //We do MHD only for the gaseous component
    
    FARGO_SAFE(UpdateMagneticField(dt,1,0,0));
    FARGO_SAFE(UpdateMagneticField(dt,0,1,0));
    FARGO_SAFE(UpdateMagneticField(dt,0,0,1));

#if !defined(STANDARD)
    FARGO_SAFE(MHD_fargo (dt)); // Perform additional field update with uniform velocity
#endif

  } 
#endif //END MHD ---------------------------------------------------------------
}

void Transport(real dt) {

  //NOTE: V_temp IS USED IN TRANSPORT
  FARGO_SAFE(copy_velocities(V2VTEMP));
#ifdef X
#ifndef STANDARD
  FARGO_SAFE(ComputeVmed(Vx_temp)); 
#endif
#endif

  transport(dt);
  
  GiveSpecificTime (t_Hydro);
  
  if (ForwardOneStep == YES) prs_exit(EXIT_SUCCESS);
  
#ifdef MHD
  if(Fluidtype == GAS) {   // We do MHD only for the gaseous component
   *(Emfx->owner) = Emfx;  // EMFs claim ownership of their storage area
   *(Emfy->owner) = Emfy;
   *(Emfz->owner) = Emfz;
 }
#endif

}

void AlgoGas1(real dt) {

  SetupHook1(); //Setup specific hook. Defaults to empty function.
  
  //Equations of state-----------------------------------------------------------
  if (Fluidtype == GAS){
#ifdef ADIABATIC
    FARGO_SAFE(ComputePressureFieldAd());
#endif
#ifdef ISOTHERMAL
    FARGO_SAFE(ComputePressureFieldIso());
#endif
#ifdef POLYTROPIC
    FARGO_SAFE(ComputePressureFieldPoly());
#endif
  }
  if(Fluidtype==DUST) FARGO_SAFE(Reset_field(Pressure));
  //-----------------------------------------------------------------------------

#ifdef X
#ifndef STANDARD
    FARGO_SAFE(ComputeVmed(Vx)); // FARGO algorithm
#endif
#endif

  InitSpecificTime (&t_Hydro, "Eulerian Hydro (no transport) algorithms");
  
  // REGARDLESS OF WHETHER WE USE FARGO, Vx IS ALWAYS THE TOTAL VELOCITY IN X
int i;

#ifdef POTENTIAL
  if(FluidIndex==0){//The potential is computed one time per timestep (FluidIndex is the local index of a fluid)
    FARGO_SAFE(compute_potential(dt));
  }
  else{
    //Correct Vtheta for a given fluid in all patches.
    if (Corotating && Current_Level == LevMax)FARGO_for_all_patches(_CorrectVtheta);
  }
#endif

#if ((defined(SHEARINGSHEET2D) || defined(SHEARINGBOX3D)) && !defined(SHEARINGBC))
  FARGO_SAFE(NonReflectingBC(Vy));
#endif
#ifdef THDIFFUSION
  if(Fluidtype==GAS){
    for(i=0;i<NSUBCYC;i++) {
      FARGO_SAFE(SubStep4_a(dt/NSUBCYC));
      FARGO_SAFE(SubStep4_b(dt/NSUBCYC));
    }
  }
#endif
//SubStep1 is divided in two halfs
//the first one performed here before,
//the second one performed afer transport
#ifdef X
  FARGO_SAFE(SubStep1_x(.5*dt));
#endif    
#ifdef Y
  FARGO_SAFE(SubStep1_y(.5*dt));
#endif  
#ifdef Z
  FARGO_SAFE(SubStep1_z(.5*dt));
#endif

#if (defined(VISCOSITY) || defined(ALPHAVISCOSITY))
  if (Fluidtype == GAS) viscosity(dt);
#endif
  
#ifndef NOSUBSTEP2
  FARGO_SAFE(SubStep2_a(dt));
  FARGO_SAFE(SubStep2_b(dt));
#endif

  // NOW: Vx INITIAL X VELOCITY, Vx_temp UPDATED X VELOCITY FROM SOURCE TERMS + ARTIFICIAL VISCOSITY

#ifdef ADIABATIC
  if(Fluidtype==GAS) FARGO_SAFE(SubStep3(dt));
#endif
  
  if(Fluidtype==DUST)FARGO_SAFE(Reset_field(Energy));
  
  GiveSpecificTime (t_Hydro);

#ifdef MHD //-------------------------------------------------------------------
  if(Fluidtype == GAS){
    InitSpecificTime (&t_Mhd, "MHD algorithms");
    FARGO_SAFE(copy_velocities(VTEMP2V));
#ifndef STANDARD // WE USE THE FARGO ALGORITHM
    FARGO_SAFE(ComputeVmed(Vx));
    FARGO_SAFE(ChangeFrame(-1, Vx, VxMed)); //Vx becomes the residual velocity
    VxIsResidual = YES;
#endif
     
    ComputeMHD(dt);

#ifndef STANDARD
    FARGO_SAFE(ChangeFrame(+1, Vx, VxMed)); //Vx becomes the total, updated velocity
    VxIsResidual = NO;
#endif 
//STANDARD
    FARGO_SAFE(copy_velocities(V2VTEMP));
// THIS COPIES Vx INTO Vx_temp
    GiveSpecificTime (t_Mhd);
  }

#endif//END MHD----------------------------------------------------------------

  InitSpecificTime (&t_Hydro, "Transport algorithms");

#if ((defined(SHEARINGSHEET2D) || defined(SHEARINGBOX3D)) && !defined(SHEARINGBC))
  FARGO_SAFE(NonReflectingBC (Vy_temp));
#endif
  
  FARGO_SAFE(copy_velocities(VTEMP2V));
}

////////////////////////////////////////////

void AlgoGas2 (real dt) {
 int i,j,k,ll;
  //At this stage V and Vtemp have the same solution,
  //except for the ghost cells, which have been properly
  //updated for V. So we copy back to Vtemp
  FARGO_SAFE(copy_velocities(V2VTEMP));
#ifdef MHD
    FARGO_SAFE(UpdateMagneticField(dt,1,0,0));
    FARGO_SAFE(UpdateMagneticField(dt,0,1,0));
    FARGO_SAFE(UpdateMagneticField(dt,0,0,1));
#endif

#if defined (MHD) && (!defined(STANDARD))
    FARGO_SAFE(MHD_fargo (dt)); // Perform additional field update with uniform velocity
#endif

#ifdef X
#ifndef STANDARD
    FARGO_SAFE(ComputeVmed(Vx_temp)); 
#endif
#endif
    
    transport(dt);

    if (Fluidtype == GAS){
#ifdef ADIABATIC
      FARGO_SAFE(ComputePressureFieldAd());
#endif
#ifdef ISOTHERMAL
      FARGO_SAFE(ComputePressureFieldIso());
#endif
#ifdef POLYTROPIC
      FARGO_SAFE(ComputePressureFieldPoly());
#endif
    }
    if(Fluidtype==DUST) FARGO_SAFE(Reset_field(Pressure));
    //Peform second half of SubStep 1
#ifdef X
    FARGO_SAFE(SubStep1_x(.5*dt));
#endif
#ifdef Y
    FARGO_SAFE(SubStep1_y(.5*dt));
#endif
#ifdef Z
    FARGO_SAFE(SubStep1_z(.5*dt));
#endif
    //Vtemp has now been updated by transport and substep1
    //We now copy this final update to V
    FARGO_SAFE(copy_velocities(VTEMP2V));

    #ifdef ACCRETION
    //For the time being let's just allow the finest level to accrete, considering that
    //it should cointaing up to the Hill sphere in it.
    if(Fluidtype==DUST && Current_Level==LevMax && ThereArePlanets)compute_accretion(dt);
    #endif
    GiveSpecificTime (t_Hydro);
  
    if (ForwardOneStep == YES) prs_exit(EXIT_SUCCESS);
  
#ifdef MHD
  if(Fluidtype == GAS) {   // We do MHD only for the gaseous component
    *(Emfx->owner) = Emfx;  // EMFs claim ownership of their storage area
    *(Emfy->owner) = Emfy;
    *(Emfz->owner) = Emfz;
  }
#endif

}
