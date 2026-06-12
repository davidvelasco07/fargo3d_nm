#include "fargo3d.h"

#define EPSILON 3.0e-6

void SolveMatrix () {
  int k;
  unsigned int count = 0;
  real omega=1.0, rho, ryx, ryz;
  real res, res0, gres, gres0;
  NormRHS ();
  res0 = reduction_full_SUM (Resid, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
#ifdef FLOAT
  MPI_Allreduce (&res0, &gres0, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Allreduce (&res0, &gres0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  gres = gres0;
  ryx = pow((Ymin(1)-Ymin(0))/(Dx*YMIN),2.0);
  ryz = pow((Ymin(1)-Ymin(0))/(Zmin(1)-Zmin(0))/YMIN,2.0);

  //Spectral radius estimate
  rho = cos(M_PI/((real)NY))+ryx*cos(M_PI/((real)NX))+ryz*cos(M_PI/((real)NZ));
  rho /= (1.0+ryx+ryz);

  //Iteration loop starts here
  while ((gres > EPSILON*gres0) && (count++ < 1000)) {
    ExecCommSame (Current_Level, ENERGYRAD);
    FARGO_SAFE (MatrixIterate (0, omega)); // First argument is parity (checkerboard)
    ExecCommSame (Current_Level, ENERGYRAD);
    FARGO_SAFE (MatrixIterate (1, omega));
    res =  reduction_full_SUM (Resid, NGHY, Ny+NGHY, NGHZ, Nz+NGHZ);
#ifdef FLOAT
    MPI_Allreduce (&res, &gres, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else
    MPI_Allreduce (&res, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
//    INSPECT_REAL (gres);
//    if (isnan(gres)) {
//      Write (nbdump++);
//      FillGhosts (StandardFields () | ENERGYRAD);
//      Write (nbdump++);
//      exit (1);
//    }
    if (count==1) 
      omega=1.0/(1.0-0.5*rho*rho);
    else 
      omega=1.0/(1.0-0.25*rho*rho*omega);
  }
  ExecCommSame (Current_Level, ENERGYRAD);
  //  INSPECT_INT (count);
  
  //  Fill_GhostsX_erad ();
}
