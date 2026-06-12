#include "fargo3d.h"


void myrestart(Field *field, int n) {
  int i,j,k;
  real *f;
  char *name;
  char filename[200];
  FILE *fi;
  int origin;

  int temp;

  f = field->field_cpu;
  name = field->name;

  if(Restart == YES) {
    sprintf(filename, "%s%s%d_%d.dat", RESTARTDIR, name, n, CPU_Rank);
    fi = fopen(filename, "r");
    if(fi == NULL) {
      printf("Error reading %s\n", filename);
      exit(1);
    }
    printf("Reading %s\n", filename);
    
    for (k=NGHZ; k<Nz+NGHZ; k++) {
      for (j=NGHY; j<Ny+NGHY; j++) {
	temp = fread(f+j*(Nx+2*NGHX)+k*Stride+NGHX, sizeof(real), Nx, fi);
      }
    }
    printf("%s OK\n", filename);
    fclose(fi);
    if(Restart_Full == YES) {
      printf("Only one restart option must be enabled.\n");
      MPI_Finalize();
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(Restart_Full == YES) {
    sprintf(filename, "%s%s%d.dat", RESTARTDIR, name, n);
    fi = fopen(filename, "r");
    if(fi == NULL) {
      printf("Error reading %s\n", filename);
      exit(1);
    }
    printf("Reading %s\n", filename);
    
    origin = (z0cell)*NX*NY + (y0cell)*NX; //z0cell and y0cell are global variables.
    for (k=NGHZ; k<Nz+NGHZ; k++) {
      fseek(fi, (origin+(k-NGHZ)*NX*NY)*sizeof(real), SEEK_SET); // critical part
      for (j=NGHY; j < Ny+NGHY; j++)
	temp = fread(f+k*Stride+j*(Nx+2*NGHX)+NGHX, sizeof(real), Nx, fi);
    }
    printf("%s OK\n", filename);
    fclose(fi);
    if(Restart == YES) {
      printf("Only one restart option must be enabled.\n");
      MPI_Finalize();
    }
  }
}

void PostRestartHook() {
  
  srand48(RSEED);

  int i,j,k;

  real r,denslocal;
  real* rho = Density->field_cpu;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
  real* vz = Vz->field_cpu;
  real* bx = Bx->field_cpu;
  real* by = By->field_cpu;
  real* bz = Bz->field_cpu;
  real* cs = Energy->field_cpu;


  printf("\nInside PostRestartHook Executing POSTRESTARHOOK SEED NUMBER=%d\n",RSEED);
  printf("====================================================================\n\n");

  myrestart(Density,0);
  myrestart(Vx,0);
  myrestart(Vy,0);
  myrestart(Vz,0);
  myrestart(Bx,0);
  myrestart(By,0);
  myrestart(Bz,0);
  myrestart(Energy,0);


//  for (k=0;k<Nz+2*NGHZ;k++) {
//    for (j=0; j<Ny+2*NGHY;j++) {
//      for (i=0; i<Nx+2*NGHX; i++) {
//	r = Ymed(j);
//	denslocal = SIGMA0/(ZMAX-ZMIN)*pow(r/R0,-SIGMASLOPE);
//	vy[l] += 0.3*(drand48()-.5)*cs[l];
//	vz[l] += 0.3*(drand48()-.5)*cs[l];
//	//rho[l] += (drand48()-.5)*denslocal;//rho[l];
//      }
//    }
//  }
}
