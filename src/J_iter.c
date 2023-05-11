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
        FARGO_SAFE(Reset_field(Total_Density)); 
        MULTIFLUID(ComputeTotalDensity()); 
        //WriteField(Total_Density,level);
#ifndef STANDARD
	      if (level == 0)
	        MULTIFLUID(ComputeVmed(Vx)); // FARGO algorithm
#endif
        #ifdef THDIFFUSION
          SelectFluid(0);//For the gas
          Compute_ThermalDiffusion();
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

      #if (defined(DRAGFORCE) && defined(DUSTSIZE))
      FARGO_SAFE(DragForce_Comm());
      #endif

      #if (defined(DUSTDIFFUSION) && !defined(DUSTSIZE))
      FARGO_SAFE(DustDiffusion_Comm());
      #endif
      
      MULTIFLUID(AlgoGas1 (dt));

      #ifdef DRAGFORCE
      FARGO_SAFE(DragForce_b(dt));
      #endif

      #ifdef DUSTDIFFUSION
      FARGO_SAFE(DustDiffusion_Main(dt));
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
      #ifdef PLANET_HEATING
      if(HEATING){
        MULTIFLUID(if(Fluidtype==GAS)compute_planetheating(dt));
      }
      #endif
      #ifdef STOCKHOLM
      if(Current_Level==0){
        MULTIFLUID(StockholmBoundary(dt));
      }
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
    //This function iterates over all fluids, 
    //therefore it is not necessary to call it with MULTIFLUID
    FARGO_SAFE(ResetFluxesLevel (level+1));
    RecursiveIteration (dt/(real)TimeStepRatio[level], level+1);
    if (TimeStepRatio[level] > 1)
      RecursiveIteration (dt/(real)TimeStepRatio[level], level+1);
    
    //This allows "level" to overwrite fluxes during transport
    time_var = LevelDate[level+1]-LevelDate[level];
    ItereLevel (time_var, level);
    LevelDate[level] = LevelDate[level+1];
    //Finer data for "level" at same date
    MULTIFLUID(ExecCommDownMean (level+1,StandardFields()));
    //We now synchronize the cpu grids of this level,
    //prior to sending information to the upper level
    TrueBC(level, StandardFields());//MULTIFLUID is called inside this function
    //Coarser data for "level+1" at same date
    MULTIFLUID(ExecCommUp (level,StandardFields()));
    return time_var;
  } 
  else {/* Finest level */
    dt_cfl_loc = CourantLimitGlobal () / (real)BaseStepRatio[level];
    MPI_Allreduce (&dt_cfl_loc, &dt_cfl, 1, MPI_DOUBLE, MPI_MIN, DomainComm); // This line may be redundant (dt_cfl_loc -> dt_cfl)
    dt = (dt_cfl < dt ? dt_cfl : dt);
    ItereLevel (dt, level);
    LevelDate[level] += dt;
    return dt;
  }
}


