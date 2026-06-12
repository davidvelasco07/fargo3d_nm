/** \file main.c

Main file of the distribution. 
Manages the call to initialization
functions, then the main loop.

*/
#include "fargo3d.h"

int   begin_i = 0, NbRestart = 0;
int   InnerOutputCounter=0, StillWriteOneOutput;
real dt;
real dtemp = 0.0;

int main(int argc, char *argv[]) {
  
  int   nn,i, OutputNumber = 0;
  char  ParameterFile[MAXLINELENGTH];
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], \
		  "-sSecCDVBndovtpfamMzib") \
	  != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], 't'))
	TimeInfo = YES;
      if (strchr (argv[i], 'f'))
	ForwardOneStep = YES;
      if (strchr (argv[i], 'C')) {
	EverythingOnCPU = YES;
#ifdef GPU
	mastererr ("WARNING: Forcing execution of all functions on CPU\n");
#else
	mastererr ("WARNING: Flag -C meaningless for a CPU built\n");
#endif
      }
      if (strchr (argv[i], 'm')) {
	Merge = YES;
#ifdef WRITEGHOSTS
	mastererr ("Cannot merge outputs when dumping ghost values.\n");
	mastererr ("'make nofulldebug' could fix this problem.\n");
	mastererr ("Starting without the -m flag could be another solution.\n");
	prs_exit (1);
#endif
      }
      if (strchr (argv[i], 's')) {
	Restart = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'p')) {
	PostRestart = YES;
	i++;
      }
      if (strchr (argv[i], 'S')) {
	Restart_Full = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'V')) {
	Dat2vtk = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect file number\n");
	  PrintUsage (argv[0]);	  
	}
      }
      if (strchr (argv[i], 'B')) {
	Vtk2dat = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect file number\n");
	  PrintUsage (argv[0]);	  
	}
      }
      if (strchr (argv[i], 'D')) {
	i++;
	DeviceManualSelection = atoi(argv[i]);
      }
    }
    else strcpy (ParameterFile, argv[i]);
  }
  
  if (ParameterFile[0] == 0) PrintUsage (argv[0]);

#ifdef MPICUDA
  EarlyDeviceSelection();
#endif
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);


#ifndef MPICUDA
  SelectDevice(CPU_Rank);
#endif
  InitVariables ();
  MPI_Barrier(MPI_COMM_WORLD);
  ReadDefaultOut ();
  ReadVarFile (ParameterFile);
  SubsDef (OUTPUTDIR, DefaultOut);
  MakeDir(OUTPUTDIR); /*Creates the output directory*/

#if !defined(X)
  NX = 1;
#endif
#if !defined(Y)
  NY = 1;
#endif
#if !defined(Z)
  NZ = 1;
#endif


  SelectWriteMethod();

#if !defined(Y) && !defined(Z)
  if (CPU_Rank==1){
    prs_error ("You cannot split your 1D mesh in x. Sequential runs only!");
  }
  if (CPU_Number > 1) {
    MPI_Finalize();
    prs_exit(EXIT_FAILURE);
  }
#endif
    
  ListVariables ("variables.par"); //Writes variables readed
  ListVariablesIDL ("IDL.var");
  ChangeArch(); /*Changes the name of the main functions
		  ChangeArch adds _cpu or _gpu if GPU is activated.*/
  split(&Gridd); /*Spliting mesh*/
  InitSpace();
  WriteDim();
  InitSurfaces();
  LightGlobalDev(); /*Filling the device global memory with light arrays.*/
  CreateFields(); // Allocate and init all fields.

  Sys = InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();
  if(Corotating)
    OMEGAFRAME = GetPsysInfo(FREQUENCY);

  CondInit();

  if(Restart == YES || Restart_Full == YES) {
    begin_i = RestartSimulation(NbRestart);
    if(strcmp(PLANETCONFIG,"NONE") != 0) {
      PhysicalTime  = GetfromPlanetFile (NbRestart, 9, 0);
      OMEGAFRAME  = GetfromPlanetFile (NbRestart, 10, 0);
      RestartPlanetarySystem (NbRestart, Sys);
    }
  }
  
  FillGhosts(PrimitiveVariables()); //Entering to the loop with all filled.
#ifdef STOCKHOLM 
    FARGO_SAFE(init_stockholm());
#endif
  for (i = begin_i; i<=NTOT; i++) { // NTOT LOOP
    if (NINTERM * (TimeStep = (i / NINTERM)) == i) {
#ifdef MHD
      FARGO_SAFE(ComputeDivergence(Bx, By, Bz));
#endif

      WriteOutputsAndDisplay();
      if (strcmp (PLANETCONFIG, "NONE") != 0)
	WritePlanetSystemFile(TimeStep);
#ifdef FARGO_COLOR
      if(CPU_Master) printf("\033[1;32m OUTPUTS %d at Physical Time t = %f OK\n\033[0m", TimeStep, PhysicalTime);
#else
      if(CPU_Master) printf("OUTPUTS %d at Physical Time t = %f OK\n", TimeStep, PhysicalTime);
#endif
      
      //      Monitor (MASS | FLUX_X | MEANBX | MEANBZ);
      ComputeMass();
      //      ComputeMomenta();
      if (TimeInfo == YES)
	GiveTimeInfo (TimeStep);
    }
    //        WritePlanetSystemFile(TimeStep);
    AlgoGas(&PhysicalTime);

    Explode();

    MonitorGlobal (MONITOR2D      | MONITORY | MONITORY_RAW|\
		   MONITORSCALAR  | MONITORZ | MONITORZ_RAW);
  }

  MPI_Finalize();  
  printf("End of simulation!\n");
  return 0;  
}
