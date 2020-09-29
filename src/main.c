/** \file main.c

Main file of the distribution. 
Manages the call to initialization
functions, then the main loop.

*/
#include "fargo3d.h"

int   begin_i = 0;
int   InnerOutputCounter=0, StillWriteOneOutput;
real dt;
real dtemp = 0.0;

int main(int argc, char *argv[]) {
  
  int   i=0, OutputNumber = 0, d;
  char  sepline[]="===========================";
  tGrid_CPU *grid;
  FluidPatch *fluid;
  sprintf (FirstCommand, "%s", argv[0]);
  sprintf (CommandLine, "%s ", argv[0]);
  while (++i < argc) {
    strncat (CommandLine, argv[i], 1023);
    strncat (CommandLine, " ", 1023);
  }

#ifdef LONGSUMMARY
  strncpy (StickyOptions, ExtractFromExecutable (YES, "", 1), 1023);
  strncpy (BoundaryFile, ExtractFromExecutable (YES, "", 3), 4095);
#endif

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
		  "-tofCmkspSVBD0#") \
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

  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_World_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_World_Number);

  NCPUS_per_domain = CPU_World_Number/NDOMAINS;
  NFluids_per_rank = NFLUIDS/NCPUS_per_domain;
  
  if (CPU_World_Number%NDOMAINS != 0) {
    printf("Error!!! The total number of Ranks (%d) must be divisible by the number of domains (%d).", CPU_World_Number, NDOMAINS);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (CPU_World_Number > NDOMAINS*NFLUIDS) {
    printf("Error!!! The total number of Ranks (%d) can not exeed NDOMAINSxNFLUIDS (%d).", CPU_World_Number, NDOMAINS*NFLUIDS);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  
  DomainColor = CPU_World_Rank%NDOMAINS;
  FluidColor  = CPU_World_Rank/NDOMAINS;
  
  MPI_Comm_split(MPI_COMM_WORLD, FluidColor, CPU_World_Rank, &DomainComm);
  MPI_Comm_rank(DomainComm, &CPU_RowRank  );
  MPI_Comm_size(DomainComm, &CPU_RowNumber);
  
  
  MPI_Comm_split(MPI_COMM_WORLD, DomainColor, CPU_World_Rank, &FluidsComm);
  MPI_Comm_rank(FluidsComm, &CPU_ColumnRank  );
  MPI_Comm_size(FluidsComm, &CPU_ColumnNumber);
  
  printf("CPU_Rank=%d\tDomainColor=%d\tFluidColor =%d\tNFLUIDS_PER_RANK=%d\n",CPU_World_Rank, DomainColor, FluidColor , NFluids_per_rank);
  
  CPU_Rank   = CPU_RowRank;
  CPU_Number = CPU_RowNumber;
  
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);

  if (strlen(xstr(VERSION)) < 2)
    sprintf (VersionString, "FARGO3D Public version 1.3");
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
  ListVariablesIDL ("IDL.var");
  ChangeArch(); /*Changes the name of the main functions
		  ChangeArch adds _cpu or _gpu if GPU is activated.*/
  
  Adapt_for_JUPITER (ParameterFile);
  
  //split(&Gridd); /*Split mesh over PEs*/
  //InitSpace();
  //WriteDim();
  //InitSurfaces();
  //LightGlobalDev(); /* Copy light arrays to the device global memory */
  //CreateFields(); // Allocate all fields.

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
    CondInit (); //Needed even for restarts: some setups have custom
		 //definitions (eg potential for setup MRI) or custom
		 //scaling laws (eg. setup planetesimalsRT).

    MULTIFLUID( begin_i  = RestartSimulation(NbRestart));
    
    if (ThereArePlanets) {
      PhysicalTime  = GetfromPlanetFile (NbRestart, 9, 0);
      OMEGAFRAME  = GetfromPlanetFile (NbRestart, 10, 0);
      RestartPlanetarySystem (NbRestart, Sys);
    }
  }
  else {
    if (ThereArePlanets)
      EmptyPlanetSystemFiles ();
      NESTEDMESHES(CondInit(););
    //CondInit(); // Initialize set up
    // Note: CondInit () must be called only ONCE (otherwise some
    // custom scaling laws may be applied several times).
  }
 
  for(i=0;i<Ngrids;i++){
    Grid_item = Grid_CPU_list; 
    do {   
      if (Grid_item->cpu == CPU_Rank && Grid_item->parent == i) { 
        printf("Grid %d\n",Grid_item->parent);
        AdaptFieldsFromJ (Grid_item);
        SelectFluid(0);
       __WriteField(Density,0);
      } 
      Grid_item= Grid_item->next;
    } while (Grid_item != NULL); 
  }
  
  exit(0);

  if (StretchOldOutput == YES) {
    StretchOutput (StretchNumber);
  }
  
  MULTIFLUID(comm(ENERGY)); //Very important for isothermal cases!

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

  GetHostsList ();
  DumpToFargo3drc(argc, argv);

  SelectArchFileName ();
#ifdef LONGSUMMARY
  ExtractFromExecutable (NO, ArchFile, 2);
#endif

  
  MULTIFLUID(FillGhosts(PrimitiveVariables()));

#ifdef STOCKHOLM 
  FARGO_SAFE(init_stockholm()); //ALREADY IMPLEMENTED MULTIFLUID COMPATIBILITY
#endif
  
#ifdef GHOSTSX
  masterprint ("\n\nNew version with ghost zones in X activated\n");
#else
  masterprint ("Standard version with no ghost zones in X\n");
#endif

  for (i = begin_i; i<=NTOT; i++) { // MAIN LOOP
    if (NINTERM * (TimeStep = (i / NINTERM)) == i) {

#if defined(MHD) && defined(DEBUG)
      FARGO_SAFE(ComputeDivergence(Bx, By, Bz));
#endif
      if (ThereArePlanets)
	WritePlanetSystemFile(TimeStep, NO);

      
#ifndef NOOUTPUTS
      MULTIFLUID(WriteOutputs(ALL));
      
#ifdef MATPLOTLIB
      Display();
#endif

      if (FluidColor == 0) {
	if(CPU_Master)
	  printf("OUTPUTS %d at date t = %f OK\n", TimeStep, PhysicalTime);
      }
#endif
      
      if (TimeInfo == YES) GiveTimeInfo (TimeStep);

    }
    
    if (NSNAP != 0) {
      if (NSNAP * (TimeStep = (i / NSNAP)) == i) {
	MULTIFLUID(WriteOutputs(SPECIFIC));
#ifdef MATPLOTLIB
	Display();
#endif
      }
    }

    if (i==NTOT)
      break;
    
    dtemp = 0.0;
    
    while (dtemp<DT) { // DT LOOP
      
      /// AT THIS STAGE Vx IS THE INITIAL TOTAL VELOCITY IN X
#ifdef X
#ifndef STANDARD
      MULTIFLUID(ComputeVmed(Vx)); // FARGO algorithm -- very important to have it here!
#endif
#endif
      /// NOW THE 2D MESH VxMed CONTAINS THE AZIMUTHAL AVERAGE OF Vx in X

#ifdef FLOOR
      MULTIFLUID(Floor());
#endif

#ifdef MHD
#ifdef OHMICDIFFUSION
      FARGO_SAFE(OhmicDiffusion_coeff());
#endif
#ifdef AMBIPOLARDIFFUSION
      FARGO_SAFE(AmbipolarDiffusion_coeff());
#endif
#ifdef HALLEFFECT
      FARGO_SAFE(HallEffect_coeff());
#endif
#endif

      // CFL condition is applied below ----------------------------------------
      MULTIFLUID(cfl());
      
      CflFluidsMin(); /*Fills StepTime with the " global min " of the
			cfl, computed from each fluid.*/
      dt = StepTime; //cfl works with the 'StepTime' global variable.
      
      dtemp+=dt;
      if(dtemp>DT)  dt = DT - (dtemp-dt); //updating dt
      //------------------------------------------------------------------------
      
      //------------------------------------------------------------------------
      /* We now compute the total density of the mesh. We need first
	 reset an array and then fill it by adding the density of each
	 fluid */
#ifdef POTENTIAL
      FARGO_SAFE(Reset_field(Total_Density)); 
      MULTIFLUID(ComputeTotalDensity());      
      MPI_Iallreduce(MPI_IN_PLACE, Total_Density->field_cpu, Nx*(Ny+2*NGHY)*(Nz+2*NGHZ),
		     MPI_DOUBLE, MPI_SUM, FluidsComm, &RequestTotalDensity);
#endif
      //------------------------------------------------------------------------
      
#ifdef COLLISIONPREDICTOR
      FARGO_SAFE(Collisions(0.5*dt, 0)); // 0 --> V is used and we update v_half.
#endif
      
      //MULTIFLUID(Sources(dt)); //v_half is used in the R.H.S

#ifdef DRAGFORCE
      FARGO_SAFE(DragForce(dt));
#endif
#ifdef DRAGFORCEALL
      FARGO_SAFE(Collisions(dt, 1));
#endif
      
#ifdef DUSTDIFFUSION
      FARGO_SAFE(DustDiffusion_Main(dt));
#endif
      
      //MULTIFLUID(Transport(dt));

      PhysicalTime+=dt;
      Timestepcount++;

#ifdef STOCKHOLM
      MULTIFLUID(StockholmBoundary(dt));
#endif

      //We apply comms and boundaries at the end of the step
      MULTIFLUID(FillGhosts(PrimitiveVariables()));


      if (FluidColor == 0) {
	if(CPU_Master) {
	  if (FullArrayComms)
	    printf("%s", "!");
	  else {
	    if (ContourComms)
	      printf("%s", ":");
	    else
	      printf("%s", ".");
	  }
	}
      }
#ifndef NOFLUSH
      fflush(stdout);
#endif
      FullArrayComms = 0;
      ContourComms = 0;
    }

    if (FluidColor == 0) {
      if(CPU_Master)
	printf("%s", "\n");
    }
    
    MULTIFLUID(MonitorGlobal (MONITOR2D      |	\
			      MONITORY       |	\
			      MONITORY_RAW   |	\
			      MONITORSCALAR  |	\
			      MONITORZ       |	\
			      MONITORZ_RAW));

    if (ThereArePlanets) {
      WritePlanetSystemFile(TimeStep, YES);
      SolveOrbits (Sys);
    }
  }
  
  MPI_Finalize();
  
  masterprint("End of the simulation!\n");
  return 0;  
}
