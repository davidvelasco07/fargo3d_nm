#include "fargo3d.h"

void AllocZmean (Field *F) {
  double *ptr;
  ptr = (double *)malloc(Nx*Ny*sizeof(double));
  if (ptr == NULL) {
    mastererr ("Out of memory in %s at line %d", __FILE__, __LINE__);
    exit (1);
  }
  F->zmean = ptr;
  ResetZmean (F);
}

int ResetZmean (Field *F) {
  int ll;
  if (F->zmean == NULL)
    AllocZmean (F);
  for (ll = 0; ll < Ny*Nx; ll++)
    F->zmean[ll] = 0.0;
  F->zp.ResetDate = PhysicalTime;
  F->zp.PreviousDate = PhysicalTime;
  //  printf ("resetting at %f\n", PhysicalTime);
}

void OutputZmean (Field *F, int counter, boolean timeaverage) {
  char filename[MAXLINELENGTH];
  FILE *out;
  int ll;
  real factor;
  if (timeaverage == YES) {
    factor = 1./(F->zp.PreviousDate-F->zp.ResetDate);
    for (ll = 0; ll < Ny*Nx; ll++)
      F->zmean[ll] *= factor;
  }
  sprintf (filename, "%s%s%d_zmean.dat", OUTPUTDIR, F->name, counter);
  masterprint ("Dumping %s%d_zmean.dat\n", F->name, counter);
  out = fopen (filename, "w");
  fwrite (F->zmean, sizeof(double), Nx*Ny, out);
  fclose (out);
  if (timeaverage == YES)
    printf ("Writing Output averaged over %f ==> %f\n", F->zp.ResetDate, F->zp.PreviousDate);
}

void Zmean (Field *F, boolean timeaverage)  {
  int i,j,k,mxy;
  real *f, t;
  real factor=1.0;
  f = F->field_cpu;
  if (F->zmean == NULL)
    AllocZmean (F);
  INPUT (F);
  t = F->zp.PreviousDate;
  if (timeaverage == YES)
    factor = PhysicalTime - t;
  for (k = NGHZ; k < Nz+NGHZ; k++) {//Vertical sum excluding ghosts
    for (j = NGHY; j < Ny+NGHY; j++) {
      for (i = 0; i < Nx; i++) {
	mxy = i+(j-NGHY)*Nx;
	F->zmean[mxy] += (double)f[l]*factor;
      }
    }
  }
  F->zp.PreviousDate = PhysicalTime;
  if (timeaverage == YES) printf ("Integrated over interval %f --> %f\n", t, PhysicalTime);
} 

void SnapshotZmean (Field *F, int counter) {
  ResetZmean (F);
  Zmean (F, NO);
  OutputZmean (F, counter, NO);
}

void DumpSnapshots (int nn) {
    SnapshotZmean (Vx, nn);
    SnapshotZmean (Vy, nn);
    SnapshotZmean (Vz, nn);
    SnapshotZmean (Bx, nn);
    SnapshotZmean (By, nn);
    SnapshotZmean (Bz, nn);
    SnapshotZmean (Density, nn);
}
