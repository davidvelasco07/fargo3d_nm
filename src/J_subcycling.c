#include "fargo3d.h"

void EvaluateLevelCost ()
{
  tGrid *grid = GridList;
  long i, nx, ny, nz;
  for (i = 0; i <= LevMax; i++)
    LevelCost[i] = 0;
  while (grid != NULL) {
    i = grid->level;
    nx = grid->gncell[0];
    ny = grid->gncell[1];
    nz = grid->gncell[2];
    LevelCost[i] += nx*ny*nz;
    grid = grid->next;
  }
}

void EvaluateBaseStepRatio ()
{
  long i;
  BaseStepRatio[0] = 1;
  for (i = 1; i <= LevMax; i++)
    BaseStepRatio[i] = BaseStepRatio[i-1]*\
      TimeStepRatio[i-1];
}

real SubCycleCost ()
{
  real cost=0.0;
  long i;
  real dt_min=1e20;
  for (i = 0; i <= LevMax; i++) {
    if (dt_min > BaseStepRatio[i]*MaxLowLevDT[i])
      dt_min = BaseStepRatio[i]*MaxLowLevDT[i];
  }
  for (i = 0; i <= LevMax; i++) {
    cost += BaseStepRatio[i]*LevelCost[i];
  }
  return cost/dt_min;
}

void FindBestSubcycling ()
{
  long i, j, remainder, nb_cases;
  real foo, cost, lowest_cost=1e28;
  real cost_nosubcycling=1.0, cost_fullsubcycling=1.0;
  char string[REFINEUPPERLIMIT+2], beststring[REFINEUPPERLIMIT+2];
  if (SubCycling == NONE) {
    for (i = 0; i < LevMax; i++)
      TimeStepRatio[i] = 1;
    EvaluateBaseStepRatio ();
    return;
  }
  if (SubCycling == FULL) {
    for (i = 0; i < LevMax; i++)
      TimeStepRatio[i] = 2;
    EvaluateBaseStepRatio ();
    return;
  }
  foo = CourantLimitGlobal ();
  EvaluateLevelCost ();
  nb_cases = (1 << LevMax);
  for (j = 0; j < nb_cases; j++) {
    remainder = j;
    for (i = 0; i < LevMax; i++) {
      string[i] = (remainder >= (1 << (LevMax-1-i)) ? '1' : '0');
      remainder %= (1 << (LevMax-1-i));
    }
    string[i] = 0;
    for (i = 0; i < LevMax; i++)
      TimeStepRatio[i] = (long)(string[i]-'0')+1;
    EvaluateBaseStepRatio ();
    cost = SubCycleCost ();
    if (j == 0) cost_nosubcycling = cost;
    if (j == nb_cases-1) cost_fullsubcycling = cost;
    if (cost < lowest_cost) {
      lowest_cost = cost;
      strncpy (beststring, string, REFINEUPPERLIMIT+1);
    }
  }
  for (i = 0; i < LevMax; i++)
    TimeStepRatio[i] = (long)(beststring[i]-'0')+1;
  EvaluateBaseStepRatio ();
  pInfo ("%s2\t", beststring);
  pInfo ("%g\t", cost_nosubcycling/lowest_cost);
  pInfo ("%g\n", cost_fullsubcycling/lowest_cost);
}
