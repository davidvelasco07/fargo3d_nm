#include "fargo3d.h"

void CflFluidsMin() {
  int i;
  real step = 1e30;
  real min;
  
  for (i=0;i<NFluids_per_rank;i++) {
    if (step > Min[i])
      step = Min[i];
  }
  
#ifdef FLOAT
  MPI_Allreduce(&step, &StepTime, 1, MPI_FLOAT, MPI_MIN, DomainComm);
#else
  MPI_Allreduce(&step, &StepTime, 1, MPI_DOUBLE, MPI_MIN, DomainComm);
#endif
  if(StepTime <= SMALLTIME) {
    masterprint("Error!!!--> Null dt\n");
    prs_exit(1);
  }
}
