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
  
  int   nn = 0, i, OutputNumber = 0;
  char  ParameterFile[MAXLINELENGTH];
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '+') {
      if (strspn (argv[i], \
		  "+#") \
	  != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], '#')) {
	i++;
	ArrayNb = atoi(argv[i]);
	EarlyOutputRename = YES;
	if (ArrayNb <= 0) {
	  masterprint ("Incorrect Array number after +# flag\n");
	  PrintUsage (argv[0]);
	}
      }
    }
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], \
		  "-tfCmspSVBD0#") \
	  != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], 't'))
	TimeInfo = YES;
      if (strchr (argv[i], 'f'))
	ForwardOneStep = YES;
      if (strchr (argv[i], '0'))
	OnlyInit = YES;
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
      if (strchr (argv[i], '#')) {
	i++;
	ArrayNb = atoi(argv[i]);
	if (ArrayNb <= 0) {
	  masterprint ("Incorrect Array number after -# flag\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'p')) {
	PostRestart = YES;
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
	Restart_Full = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect output number\n");
	  PrintUsage (argv[0]);	  
	}
      }
      if (strchr (argv[i], 'B')) {
	Vtk2dat = YES;
	Restart_Full = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect output number\n");
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

  /* This must be placed ***BEFORE*** reading the input files in case of a restart */
  if ((ArrayNb) && (EarlyOutputRename == YES)) {
    i = strlen(OUTPUTDIR);
    if (OUTPUTDIR[i-1] == '/') OUTPUTDIR[i-1] = 0;//Remove trailing slash if any
    sprintf (OUTPUTDIR, "%s%06d/", OUTPUTDIR, ArrayNb); //Append numerical suffix
    /* There is no need to perform the wildcard (@) substitution. This has already been done */
    printf ("\n\n***\n\nNew Output Directory is %s\n\n***\n\n", OUTPUTDIR);
    MakeDir(OUTPUTDIR); /*Create the output directory*/
  }


  MakeDir(OUTPUTDIR); /*Create the output directory*/

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
    prs_error ("You cannot split a 1D mesh in x. Sequential runs only!");
  }
  if (CPU_Number > 1) {
    MPI_Finalize();
    prs_exit(EXIT_FAILURE);
  }
#endif
    
  ListVariables ("variables.par"); //Writes all variables defined in set up
  ListVariablesIDL ("IDL.var");
  ChangeArch(); /*Changes the name of the main functions
		  ChangeArch adds _cpu or _gpu if GPU is activated.*/
  split(&Gridd); /*Split mesh over PEs*/
  InitSpace();
  WriteDim();
  InitSurfaces();
  LightGlobalDev(); /* Copy light arrays to the device global memory */
  CreateFields(); // Allocate all fields.

  Sys = InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();
  if(Corotating)
    OMEGAFRAME = GetPsysInfo(FREQUENCY);

  CondInit(); // Initialize set up

  if(Restart == YES || Restart_Full == YES) {
    begin_i = RestartSimulation(NbRestart);
    if(strcmp(PLANETCONFIG,"NONE") != 0) {
      PhysicalTime  = GetfromPlanetFile (NbRestart, 9, 0);
      PhysicalTimeRestart = PhysicalTime;
      OMEGAFRAME  = GetfromPlanetFile (NbRestart, 10, 0);
      RestartPlanetarySystem (NbRestart, Sys);
    }
  }
  /* This must be placed ***after*** reading the input files in case of a restart */
  if ((ArrayNb) && (EarlyOutputRename == NO)) {
    i = strlen(OUTPUTDIR);
    if (OUTPUTDIR[i-1] == '/') OUTPUTDIR[i-1] = 0;//Remove trailing slash if any
    sprintf (OUTPUTDIR, "%s%06d/", OUTPUTDIR, ArrayNb); //Append numerical suffix
    /* There is no need to perform the wildcard (@) substitution. This has already been done */
    printf ("\n\n***\n\nNew Output Directory is %s\n\n***\n\n", OUTPUTDIR);
    MakeDir(OUTPUTDIR); /*Create the output directory*/
    ListVariables ("variables.par"); //Writes all variables defined in set up
    ListVariablesIDL ("IDL.var");
    InitSpace();
    WriteDim ();
  }
  
  FillGhosts(PrimitiveVariables()); 
#ifdef STOCKHOLM 
  FARGO_SAFE(init_stockholm());
#endif
  DumpSnapshots (nn);
  for (i = begin_i; i<=NTOT; i++) { // MAIN LOOP
    if (NINTERM * (TimeStep = (i / NINTERM)) == i) {
#if defined(MHD) && defined(DEBUG)
      FARGO_SAFE(ComputeDivergence(Bx, By, Bz));
#endif

#ifndef NOOUTPUTS
      WriteOutputsAndDisplay();

      if (strcmp (PLANETCONFIG, "NONE") != 0)
	WritePlanetSystemFile(TimeStep);
#ifdef FARGO_COLOR
      if(CPU_Master) printf("\033[1;32m OUTPUTS %d at date t = %f OK\n\033[0m", TimeStep, PhysicalTime);
#else
      if(CPU_Master) printf("OUTPUTS %d at date t = %f OK\n", TimeStep, PhysicalTime);
#endif
#endif
      

      //      ComputeMass();
      if (TimeInfo == YES)
	GiveTimeInfo (TimeStep);
    }
    AlgoGas(&PhysicalTime);
    MonitorGlobal (MONITOR2D      | MONITORY | MONITORY_RAW|\
		   MONITORSCALAR  | MONITORZ | MONITORZ_RAW);
    nn++;
    DumpSnapshots (nn);
  }

  MPI_Finalize();  
  printf("End of simulation!\n");
  return 0;  
}
