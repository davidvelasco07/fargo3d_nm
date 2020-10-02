#include "fargo3d.h"

int StandardFields () {
  int options;
  options = DENS;
#ifdef X
  options |= VX;
#endif
#ifdef Y
  options |= VY;
#endif
#ifdef Z
  options |= VZ;
#endif
#ifdef ADIABATIC
  options |= ENERGY;
#endif
#ifdef LABELED
  options |= LABEL;
#endif
  return options;
}

void TrueBC (long lev, int options){
  int dim, side, dimp1, dimp2, k, size;
  tGrid_CPU *item;
  item = Grid_CPU_list;
  
  do {
    if (item->cpu == CPU_Rank) {
      if (item->level == lev) {
	      FARGO_SAFE(AdaptFieldsFromJ (item));
	      MULTIFLUID(FillGhosts(options)); //does a comm and boundaries() 
      }
    }
    item = item->next;
  } while (item != NULL);

}

void UpdateCourantLimit (long level)
{
  tGrid_CPU *item;
  FluidPatch *Fluid;
  real dt;
  item = Grid_CPU_list;
  MaxLowLevDT_loc[level];
  while (item != NULL) {
    if (item->cpu == CPU_Rank) {
      if (level == item->level) {
        FARGO_SAFE(AdaptFieldsFromJ (item));
#ifndef STANDARD
	      if (level == 0)
	        MULTIFLUID(ComputeVmed(Vx)); // FARGO algorithm
#endif
        MULTIFLUID(cfl());
        CflFluidsMin();
	      //cfl works with the 'StepTime' global variable.
        MaxLowLevDT_loc[level] = StepTime;
      }
    }
    item = item->next;
  } /* Note that at this stage the MaxLowLevDT depends on the p.e. */
}



real CourantLimitGlobal ()
{
  long i;
  real dt_min;
  for (i = 0; i <= LevMax; i++) {
    if (LevelHasChangedSinceCFL[i] == YES) {
      UpdateCourantLimit(i);
      LevelHasChangedSinceCFL[i] = NO;
    }
  }
  MPI_Allreduce (MaxLowLevDT_loc, MaxLowLevDT, LevMax+1, MPI_DOUBLE, MPI_MIN, DomainComm);
  dt_min = BaseStepRatio[0]*MaxLowLevDT[0];
  for (i = 1; i <= LevMax; i++) {
    if (dt_min > BaseStepRatio[i]*MaxLowLevDT[i])
      dt_min = BaseStepRatio[i]*MaxLowLevDT[i];
  }
  return (dt_min < DT ? dt_min : DT);
}
 
void ItereLevel (real dt, long level)
{
  tGrid_CPU *item;
  FluidPatch *Fluid;
  
  TrueBC(level,StandardFields());
  
  item = Grid_CPU_list;
  while (item != NULL) {//Source 1st half 
    if ((level == item->level) && (item->cpu == CPU_Rank)) {
      //To make this work Adapt... has to Change Fluids-> to point to the fluids of the grid
      //We also have to include Vi_temp as a particular field for each fluid and each patch
      FARGO_SAFE(AdaptFieldsFromJ (item));
      #ifdef POTENTIAL
      FARGO_SAFE(Reset_field(Total_Density)); 
      MULTIFLUID(ComputeTotalDensity());      
      MPI_Iallreduce(MPI_IN_PLACE, Total_Density->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		     MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestTotalDensity);
      #endif
      MULTIFLUID(AlgoGas1 (dt));
      
      #ifdef DRAGFORCE
      //FARGO_SAFE(DragForce(.5*dt));
      #endif
    }
    item = item->next;
  }
  
  TrueBC(level,StandardFields());//This may be needed with a grid of processes
  
  item = Grid_CPU_list;
  while (item != NULL) {//Transport + Source 2nd half
    if ((level == item->level) && (item->cpu == CPU_Rank)) {
      FARGO_SAFE(AdaptFieldsFromJ (item));
      MULTIFLUID(AlgoGas2 (dt));
      #ifdef DRAGFORCE
      //FARGO_SAFE(DragForce(.5*dt));
      #endif
    }
    item = item->next;
  }

  LevelHasChangedSinceCFL[level] = YES;
}

real RecursiveIteration (real dt, long level)
{
  real time_var, dt_cfl_loc, dt_cfl;
  long i;
  if (level == 0) {		/* We begin a global timestep */
    for (i = 0; i <= LevMax; i++) {
      LevelHasChangedSinceCFL[i] = YES;
    }
    FARGO_SAFE(FindBestSubcycling ());
  }
  if (level < LevMax) {/* Not finest level */
    FARGO_SAFE(ResetFluxesLevel (level+1));
    RecursiveIteration (dt/(real)TimeStepRatio[level], level+1);
    if (TimeStepRatio[level] > 1)
      RecursiveIteration (dt/(real)TimeStepRatio[level], level+1);
    
    //This allows "level" to overwrite fluxes during transport
    time_var = LevelDate[level+1]-LevelDate[level];
    ItereLevel (time_var, level);
    LevelDate[level] = LevelDate[level+1];
   
    MULTIFLUID(ExecCommDownMean (level+1,StandardFields()));

    //Finer data for "level" at same date
    //MULTIFLUID is called inside this function for FillGhost
    //Which then calls CommSame and boundaries
    TrueBC(level, StandardFields());
    //Includes comm same and boundaries 
    
    MULTIFLUID(ExecCommUp (level,StandardFields()));
    //Coarser data for "level+1" at same date

    return time_var;
  } else {/* Finest level */
    dt_cfl_loc = CourantLimitGlobal () / (real)BaseStepRatio[level];
    MPI_Allreduce (&dt_cfl_loc, &dt_cfl, 1, MPI_DOUBLE, MPI_MIN, DomainComm);
    dt = (dt_cfl < dt ? dt_cfl : dt);
    ItereLevel (dt, level);
    LevelDate[level] += dt;
    return dt;
  }
}


