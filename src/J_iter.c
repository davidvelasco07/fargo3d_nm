#include "fargo3d.h"
#include "J_jupiter.h"

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
	FillGhosts(options); //does a comm and boundaries() 
      }
    }
    item = item->next;
  } while (item != NULL);

}

void UpdateCourantLimit (level)
     long level;
{
  tGrid_CPU *item;
  FluidPatch *Fluid;
  real dt;
  item = Grid_CPU_list;
  MaxLowLevDT_loc[level] = 1e20;
  while (item != NULL) {
    if (item->cpu == CPU_Rank) {
      if (level == item->level) {
        FARGO_SAFE(AdaptFieldsFromJ (item));
#ifndef STANDARD
	if (level == 0)
	  FARGO_SAFE(ComputeVmed(Vx)); // FARGO algorithm
#endif
        FARGO_SAFE(cfl());
	
	dt = step_time; //cfl works with the 'step_time' global variable.
	//if(CPU_Rank == 0)
	//  printf ("Limite dt=%.14g sobre el nivel %ld\n", dt, level);
	if (dt < MaxLowLevDT_loc[level]) MaxLowLevDT_loc[level] = dt;
      }
    }
    item = item->next;
  } /* Note that at this stage the MaxLowLevDT depends on the p.e. */
}



real CourantLimitGlobal ()
{
  long i;
  real dt_min = 1e20;
  for (i = 0; i <= LevMax; i++) {
    if (LevelHasChangedSinceCFL[i] == YES) {
      UpdateCourantLimit(i);
      LevelHasChangedSinceCFL[i] = NO;
    }
  }
  MPI_Allreduce (MaxLowLevDT_loc, MaxLowLevDT, LevMax+1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  for (i = 0; i <= LevMax; i++) {
    if (dt_min > BaseStepRatio[i]*MaxLowLevDT[i])
      dt_min = BaseStepRatio[i]*MaxLowLevDT[i];
  }
  return (dt_min < DT ? dt_min : DT);
}
 
void ItereLevel (dt, level)
     real dt;
     long level;
{
  tGrid_CPU *item;
  FluidPatch *Fluid;
  
  TrueBC(level,StandardFields());

  item = Grid_CPU_list;
  while (item != NULL) {//Source 1st half 
    if ((level == item->level) && (item->cpu == CPU_Rank)) {
      FARGO_SAFE(AdaptFieldsFromJ (item));
      AlgoGas (dt);
    }
    item = item->next;
  }
 
  TrueBC(level,StandardFields());//This may be needed with a grid of processes
  
  item = Grid_CPU_list;
  while (item != NULL) {//Transport + Source 2nd half
    if ((level == item->level) && (item->cpu == CPU_Rank)) {
      FARGO_SAFE(AdaptFieldsFromJ (item));
      AlgoGas2 (dt);
    }
    item = item->next;
  }
  
  
#ifdef RADIATIVE_TRANSFER
  item = Grid_CPU_list;
  ExecCommUp (level-1, StandardFields() | ENERGYRAD);
  FARGO_for_all_patches_level_intarg (FillGhosts, level, StandardFields() | ENERGYRAD);
  ExecCommDownMean (level, ENERGYRAD);
  while (item != NULL) {
    if ((level == item->level) && (item->cpu == CPU_Rank)) {
      FARGO_SAFE(AdaptFieldsFromJ (item));
      FARGO_SAFE(RT_main(dt));
    }
    item = item->next;
  }
#endif
  LevelHasChangedSinceCFL[level] = YES;
}

real RecursiveIteration (dt, level)
     real dt;
     long level;
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
    FARGO_SAFE(ExecCommDownFlux (level+1));
    //This allows "level" to overwrite fluxes during transport
    time_var = LevelDate[level+1]-LevelDate[level];
    ItereLevel (time_var, level);
    LevelDate[level] = LevelDate[level+1];
   
    FARGO_SAFE(ExecCommDownMean (level+1,StandardFields()));
    //Finer data for "level" at same date
    FARGO_SAFE(TrueBC(level, StandardFields()));
    //Includes comm same and boundaries 
    
    FARGO_SAFE(ExecCommUp (level,StandardFields()));
    //Coarser data for "level+1" at same date
     
    return time_var;
  } else {/* Finest level */
    dt_cfl_loc = CourantLimitGlobal () / (real)BaseStepRatio[level];
    MPI_Allreduce (&dt_cfl_loc, &dt_cfl, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = (dt_cfl < dt ? dt_cfl : dt);
    ItereLevel (dt, level);
    LevelDate[level] += dt;
    return dt;
  }
}


