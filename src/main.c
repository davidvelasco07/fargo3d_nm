/** \file main.c

Main file of the distribution. 
Manages the call to initialization
functions, then the main loop.

This version of the code imports some features of the nested-mesh code JUPITER.

*/
#include "fargo3d.h"
#include "J_jupiter.h"

int   begin_i = 0, level;
int   InnerOutputCounter=0, StillWriteOneOutput;
real dt;
real dtemp = 0.0;

int main(int argc, char *argv[]) {
  int   i=0, OutputNumber = 0, d;
  real OMEGAFRAME_PARAMETER;
  char  sepline[]="===========================";
  real NextDate;
  tGrid_CPU *grid;
  FluidPatch *fluid;
  if (argc == 1) PrintUsage (argv[0]);

  sprintf (FirstCommand, "%s", argv[0]);
  sprintf (CommandLine, "%s ", argv[0]);
  while (++i < argc) {
    strncat (CommandLine, argv[i], 1023);
    strncat (CommandLine, " ", 1023);
  }

  //strncpy (StickyOptions, ExtractFromExecutable (YES, "", 1), 1023);
  //strncpy (BoundaryFile, ExtractFromExecutable (YES, "", 3), 1023);
  
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i+=d) {
    d=1;
    if (*(argv[i]) == '+') {
      if (strspn (argv[i], \
		  "+S#D") \
	  != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], '#')) {
	d=2;
	ArrayNb = atoi(argv[i+1]);
	EarlyOutputRename = YES;
	if (ArrayNb <= 0) {
	  masterprint ("Incorrect Array number after +# flag\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'D')) {
	d=2;
	strcpy (DeviceFile, argv[i+1]);
	DeviceFileSpecified = YES;
      }
      if (strchr (argv[i], 'S')) {
	d=2;
	StretchNumber = atoi(argv[i+1]);
	StretchOldOutput = YES;
      }
    }
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], \
		  "-tofCmMkspSVBD0#jr")		\
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
      }
      if (strchr (argv[i], 'k')) {
	Merge = NO;
      }
      if (strchr (argv[i], 'o')) {
	RedefineOptions = YES;
	ParseRedefinedOptions (argv[i+1]) ;
	d=2;
      }
      if (strchr (argv[i], 's')) {
	Restart = YES;
	d=2;
	NbRestart = atoi(argv[i+1]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'j')) {
	MaxLevel = atoi(argv[i+1]);
	if (MaxLevel < 0) masterprint ("Incorrect maximal level specification");
	d=2;
      }
      if (strchr (argv[i], 'r')) {
        AddSubPatch = YES;
	sprintf (SubPatchInfo, "%s", argv[i+1]);
	d=2;
      }
      if (strchr (argv[i], '#')) {
	d=2;
	ArrayNb = atoi(argv[i+1]);
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
	d=2;
	NbRestart = atoi(argv[i+1]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'M')) {
	d=2;
	NbRestart = atoi(argv[i+1]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect merge number\n");
	  PrintUsage (argv[0]);
	}
	MergePostProcess = YES;
      }
      if (strchr (argv[i], 'V')) {
	Dat2vtk = YES;
	Restart_Full = YES;
	d=2;
	NbRestart = atoi(argv[i+1]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect output number\n");
	  PrintUsage (argv[0]);	  
	}
      }
      if (strchr (argv[i], 'B')) {
	Vtk2dat = YES;
	Restart_Full = YES;
	d=2;
	NbRestart = atoi(argv[i+1]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect output number\n");
	  PrintUsage (argv[0]);	  
	}
      }
      if (strchr (argv[i], 'D')) {
	d=2;
	DeviceManualSelection = atoi(argv[i+1]);
      }
    }
    else strcpy (ParameterFile, argv[i]);
  }

#ifdef WRITEGHOSTS
  if (Merge == YES) {
	mastererr ("Cannot merge outputs when dumping ghost values.\n");
	mastererr ("'make nofulldebug' could fix this problem.\n");
	mastererr ("Using the -k flag could be another solution.\n");
	prs_exit (1);
  }
#endif
  

#ifdef MPICUDA
  EarlyDeviceSelection();
#endif
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);

  if (strlen(xstr(VERSION)) < 2)
    sprintf (VersionString, "FARGO3D Public version 1.2");
  else
    sprintf (VersionString, "FARGO3D git version %s", xstr(VERSION));
  masterprint("\n\n%s\n%s\nSETUP: '%s'\n%s\n\n",
	      sepline, VersionString, xstr(SETUPNAME), sepline);
  
  if ((ParameterFile[0] == 0) || (argc == 1)) PrintUsage (argv[0]);

#ifndef MPICUDA
  SelectDevice(CPU_Rank);
#endif
  InitVariables ();
  MPI_Barrier(MPI_COMM_WORLD);
  ReadDefaultOut ();
  ReadVarFile (ParameterFile);

  if (strcmp (PLANETCONFIG, "NONE") != 0)
    ThereArePlanets = YES;

  OMEGAFRAME_PARAMETER = OMEGAFRAME;

  if (ORBITALRADIUS > 1.0e-30){
    YMIN *= ORBITALRADIUS;
    YMAX *= ORBITALRADIUS;
    DT   *= sqrt(ORBITALRADIUS*ORBITALRADIUS*ORBITALRADIUS);
  }


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
  ListVariablesIDL ("IDL_params.dat");

  ChangeArch(); /*Changes the name of the main functions
		  ChangeArch adds _cpu or _gpu if GPU is activated.*/
  
  FARGO_SAFE(Adapt_for_JUPITER (ParameterFile));
  
  if (MergePostProcess)
    merge (NbRestart);
  
  //SubCycling = AUTO;
  //SubCycling = NONE;
  SubCycling = FULL;

  // LightGlobalDev(); /* Copy light arrays to the device global memory */  CHECK LATER

  Sys = InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();
  
  if(Corotating)
    OMEGAFRAME = GetPsysInfo(FREQUENCY);
  OMEGAFRAME0 = OMEGAFRAME;
  /* We need to keep track of initial azimuthal velocity to correct
the target velocity in Stockholm's damping prescription. We copy the
value above *after* rescaling, and after any initial correction to
OMEGAFRAME (which is used afterwards to build the initial Vx field. */
    
  if(Restart == YES || Restart_Full == YES) {
    //merge(NbRestart);
    FARGO_for_all_patches (CondInit); //Needed even for restarts: some setups have custom
		 //definitions (eg potential for setup MRI) or custom
		 //scaling laws (eg. setup planetesimalsRT).
    grid = Grid_CPU_list;
    while (grid != NULL){
      if (grid->cpu == CPU_Rank) {
	fluid=grid->Fluid;
	while (fluid != NULL) {
	    ReadField(fluid,NbRestart);
	    fluid = fluid->next;
	  }
      }
      grid=grid->next;
    }
    begin_i = NbRestart*NINTERM;
    //begin_i = RestartSimulation(NbRestart);
    if (ThereArePlanets) {
      PhysicalTime  = GetfromPlanetFile (NbRestart, 9, 0);
      OMEGAFRAME  = GetfromPlanetFile (NbRestart, 10, 0);
      RestartPlanetarySystem (NbRestart, Sys);
    }
  }
  else {
    if (ThereArePlanets)
      EmptyPlanetSystemFiles ();
   
    FARGO_SAFE(FARGO_for_all_patches (CondInit));
    
    // Note: CondInit () must be called only ONCE (otherwise some
    // custom scaling laws may be applied several times).
  }

  if (StretchOldOutput == YES) {
    OMEGAFRAME = OMEGAFRAME_PARAMETER;
    StretchOutput (StretchNumber);
  }

  //  FARGO_SAFE(comm(ENERGY)); //Very important for isothermal cases! CHECK LATER

  /* This must be placed ***after*** reading the input files in case of a restart */
  if ((ArrayNb) && (EarlyOutputRename == NO)) {
    i = strlen(OUTPUTDIR);
    if (OUTPUTDIR[i-1] == '/') OUTPUTDIR[i-1] = 0;//Remove trailing slash if any
    sprintf (OUTPUTDIR, "%s%06d/", OUTPUTDIR, ArrayNb); //Append numerical suffix
    /* There is no need to perform the wildcard (@) substitution. This has already been done */
    printf ("\n\n***\n\nNew Output Directory is %s\n\n***\n\n", OUTPUTDIR);
    MakeDir(OUTPUTDIR); /*Create the output directory*/
    ListVariables ("variables.par"); //Writes all variables defined in set up
    ListVariablesIDL ("IDL_params.dat");
    //    InitSpace();
    WriteDim ();
  }

  GetHostsList ();
  DumpToFargo3drc(argc, argv);

  SelectArchFileName ();
  //ExtractFromExecutable (NO, ArchFile, 2);
  
  //FillGhosts(PrimitiveVariables()); 
#ifdef STOCKHOLM 
  FARGO_for_all_patches(init_stockholm);
#endif

#ifdef GHOSTSX
  masterprint ("\n\nNew version with ghost zones in X activated\n");
#else
  masterprint ("Standard version with no ghost zones in X\n");
#endif
 
  for (level = 0; level <= LevMax; level++) {
#ifdef LUMINOUS_SRC
#if LUMINOUS_SRC == 1
    //Spherical release
    FARGO_for_all_patches_level (Count_Luminous_cells,level);
    
#endif
#endif
    FARGO_SAFE(FARGO_for_all_patches_level_intarg (FillGhosts, level,StandardFields() | ENERGY));
    FARGO_SAFE(ExecCommUp (level,StandardFields() | ENERGY));
  }
 
#ifdef STOCKHOLM 
  FARGO_for_all_patches_level (init_stockholm, 0);
#endif

  tGrid_CPU *current;
  real mass, totalmass;
  for (i = begin_i; i<=NTOT; i++) { // MAIN LOOP
    if (NINTERM * (TimeStep = (i / NINTERM)) == i) { //Every NINTERM steps we save stuff

#if defined(MHD) && defined(DEBUG)
      FARGO_SAFE(ComputeDivergence(Bx, By, Bz));
#endif
      if (ThereArePlanets)
	WritePlanetSystemFile(TimeStep, NO);
      
#ifndef NOOUTPUTS
      //WriteOutputsAndDisplay(ALL);

      FARGO_SAFE(Write (TimeStep));
   ;
#ifdef REMAP
      WriteY(TimeStep); //Dumping domain_y when the mesh is varying.
      printf("YMIN=%lf, YMAX=%lf\n",YMIN,YMAX);
#endif

      if(CPU_Master) printf("\nOUTPUTS %d at date t = %f OK\n", TimeStep, PhysicalTime);
#endif

      if (TimeInfo == YES)
	GiveTimeInfo (TimeStep);

    }//END Every NINTERM steps we save stuff
    
    if (NSNAP != 0) {
      if (NSNAP * (TimeStep = (i / NSNAP)) == i) {
	WriteOutputsAndDisplay(SPECIFIC);
      }
    }
    
    NextTime = PhysicalTime+DT;
    while (PhysicalTime < NextTime-DT/1e10){
      PhysicalTime += RecursiveIteration (NextTime-PhysicalTime, 0L);
      //exit(1);
    }
    mass=0.0;
    totalmass = 0.0;
    current = Grid_CPU_list;
    while(current != NULL){
      if (current->cpu == CPU_Rank) {
	AdaptFieldsFromJ (current);
	mass += ComputeMass(1);	
      }
      current = current->next;
    }
    MPI_Allreduce (&mass, &totalmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(CPU_Rank==0)printf("\nTotalMass = %3.17f \n", totalmass );
    if (ThereArePlanets) {
      WriteTorqueAndWork(TimeStep, 0);
      WritePlanetSystemFile(TimeStep, YES);
      SolveOrbits (Sys);
    }
    MonitorGlobal (MONITOR2D      | MONITORY | MONITORY_RAW|	\
		 MONITORSCALAR  | MONITORZ | MONITORZ_RAW);
  
    MonitorNested (MONITORSCALAR);
  }
  
  
  MPI_Finalize();  
  printf("End of simulation!\n");
  return 0;  
}
