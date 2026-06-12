#define __SWITCH_
#include "jupiter.h"
#undef __SWITCH_

char *CoordNames[] = {"X", "Y", "Z", "Azimuth", "Radius", "Z", "Azimuth", "Radius", "Co-latitude"};
char *SCoordNames[] = {"X", "Y", "Z", "Phi", "Rad", "Z", "Phi", "Rad", "Theta"};


long switches (argc, argv, parameter_file)
     int argc;
     char *argv[];
     char *parameter_file;
{
  long NbRestart = 0, i,d=1;
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (parameter_file, "");
  for (i = 1; i < argc; i+=d) {
    d=1;
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], "-acegijkmnopqrstvxzSh") != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], 'o')) {
	RedefineOptions = YES;
	ParseRedefinedOptions (argv[i+1]) ;
	d=2;
      }
      if (strchr (argv[i], 's')) {
	Restart = YES;
	NbRestart = atoi(argv[i+1]);
	if (NbRestart < 0) prs_error ("Incorrect restart number");
	d=2;
      }
      if (strchr (argv[i], 'S')) {
	Stretch = YES;
	NbRestart = atoi(argv[i+1]);
	if (NbRestart < 0) prs_error ("Incorrect restart/stretch number");
	d=2;
      }
      if (strchr (argv[i], 'j')) {
	MaxLevel = atoi(argv[i+1]);
	if (MaxLevel < 0) prs_error ("Incorrect maximal level specification");
	d=2;
      }
      if (strchr (argv[i], 'm')) {
	Merge = YES;
	NbRestart = atoi(argv[i+1]);
	if (NbRestart < 0) prs_error ("Incorrect merging number");
	d=2;
      }
      if (strchr (argv[i], 'r')) {
	AddSubPatch = YES;
	sprintf (SubPatchInfo, "%s", argv[i+1]);
	d=2;
      }
      if (strchr (argv[i], 'x')) {
	RayTracing = YES;
	Restart = YES;
	sprintf (RayTracingInfo, "%s", argv[i+1]);
	sscanf (RayTracingInfo, "%ld", &NbRestart);
	if (NbRestart < 0) prs_error ("Incorrect output number");
	d=2;
      }
      if (strchr (argv[i], 'a'))
	AllCPUs = YES;
      if (strchr (argv[i], 'p'))
	MonitorTimeSpent = YES;
      if (strchr (argv[i], 'q'))
	QuietOutput = YES;
      if (strchr (argv[i], 'n'))
	Disable = YES;
      if (strchr (argv[i], 'e'))
	EneMon = YES;
      if (strchr (argv[i], 'k'))
	TorqStock = YES;
      if (strchr (argv[i], 't'))
	TorqueMon = YES;
      if (strchr (argv[i], 'i'))
	SuperImpose = YES;
      if (strchr (argv[i], 'h'))
	SmoothTaper = YES;
      if (strchr (argv[i], 'c')) {
	CFLMonitor = YES;
	CoarseRayTracing = YES;
      }
      if (strchr (argv[i], 'z'))
	MonitorCons = NO;
      if (strchr (argv[i], 'g'))
	WriteGhosts = YES;
    }
    else strcpy (parameter_file, argv[i]);
  }
  if (parameter_file[0] == 0) PrintUsage (argv[0]);
  SubsDef ("input file", parameter_file, DefaultIn);
  ExtractPath (parameter_file, InputPath, ShortParName);
  pInfo ("Input file path is %s\n", InputPath);
  pInfo ("Parameter file name is %s\n", ShortParName);
  return NbRestart;
}

void SetGlobalVariables () {
  long i, j, bc[6];
  boolean set;
  char spacingdim[3][MAXLINELENGTH];
  char message[MAXLINELENGTH];
  FILE *TGF;
  
  /* First we check that COORDPERMUT */
  /* is indeed a permutation of 123 */
  if ((strlen(COORDPERMUT) != 3) || (strspn("123",COORDPERMUT) != 3) || (strcspn("123",COORDPERMUT))) prs_error ("COORDPERMUT should be a permutation of 123");
  if ((*(COORDTYPE+1) == 'y') || (*(COORDTYPE+1) == 'Y'))
    CoordType = CYLINDRICAL;
  if ((*(COORDTYPE+0) == 's') || (*(COORDTYPE+0) == 'S'))
    CoordType = SPHERICAL;
  for (i = 0; i < 3; i++)
    InvCoordNb[i] = (long)(*(COORDPERMUT+i))-(long)('1');
  strcpy(spacingdim[0], DIM1SPACING);
  strcpy(spacingdim[1], DIM2SPACING);
  strcpy(spacingdim[2], DIM3SPACING);
  for (i = 0; i < 3; i++) {
    if ((*(spacingdim[i]) == 'G') || (*(spacingdim[i]) == 'g'))
      SpacingDim[i] = GEOMETRIC;
  }
  /* Output directory rectification (add one trailing slash if needed) */
  if (*(OUTPUTDIR+strlen(OUTPUTDIR)-1) != '/')
    strcat (OUTPUTDIR, "/");
  pInfo ("\nThe coordinate system is ");
  switch (CoordType) {
  case CARTESIAN:   pInfo("cartesian\n"); __CARTESIAN = YES; break;
  case CYLINDRICAL: pInfo("cylindrical\n"); __CYLINDRICAL = YES; break;
  case SPHERICAL:   pInfo("spherical\n"); __SPHERICAL = YES; break;
  }
  pInfo ("\nThe coordinates are respectively : ");
  for (i = 0; i < 3; i++) {
    pInfo(CoordNames[CoordType*3+InvCoordNb[i]]);
    if (i == 1) pInfo(" and ");
    if (i == 0) pInfo(", ");
  }
  pInfo ("\n");
  for (i = 0; i < NDIM; i++) {
    sprintf (message, "The coordinate spacing along dimension %ld is ", i+1);
    switch (SpacingDim[i]) {
    case ARITHMETIC: strcat (message, "arithmetic"); break;
    case GEOMETRIC: strcat (message, "geometric"); break;
    }
    pInfo (message);
    pInfo ("\n");
  }
  /* We now invert the coord. permutation, */
  /* i.e. we take its reciprocal map */
  for (i = 0; i < 3; i++) {
    j = 0;
    while (InvCoordNb[j] != i) j++;
    CoordNb[i] = j;
  }
  _X_ = _Radial_ = _RAD_ = CoordNb[0];
  _Y_ = _Azimuthal_ = _AZIM_ = CoordNb[1];
  _Z_ = _Vertical_ = _Colatitude_ = _COLAT_ = CoordNb[2];
  if (EmbeddedGridFile == YES) {
    if (*GRIDFILE != 0) {
      pWarning ("You have defined the GRIDFILE variable");
      pWarning ("but you have also an embedded grid file, specified");
      pWarning ("via the #GRIDINFO directive in your parameter file.");
      pWarning ("I will discard the embedded information and I will try");
      pWarning ("to find the grid file, and I will issue an error");
      pWarning ("message if I cannot find it.");
    } else {
      strcpy (GRIDFILE, EmbeddedGridFileName);
    }
  }
  if (((*CONDLIM1 != 'U') || (*CONDLIM2 != 'U') || (*CONDLIM3 != 'U')) &&\
      (*GRIDFILE != 0)) {
    pWarning ("You have specified boundary conditions via CONDLIM variables,");
    pWarning ("but you have also defined a grid file. I will use the boundary conditions");
    pWarning ("specified in that grid file, and I will disregard whatever would be");
    pWarning ("imposed by the CONDLIM variables.");
  }
  if (*GRIDFILE == 0) {	/* The template is then written in a temporary grid file */
    switch (*CONDLIM1) {
    case 'O': bc[0] = bc[3] = 2; break;	/* outflow */
    case 'o': bc[0] = bc[3] = 2; break;
    case '2': bc[0] = bc[3] = 2; break;
    case 'R': bc[0] = bc[3] = 1; break;	/* reflecting */
    case 'r': bc[0] = bc[3] = 1; break;
    case '1': bc[0] = bc[3] = 1; break;
    default : bc[0] = bc[3] = 0;
    }
    switch (*CONDLIM2) {
    case 'O': bc[1] = bc[4] = 2; break;
    case 'o': bc[1] = bc[4] = 2; break;
    case '2': bc[1] = bc[4] = 2; break;
    case 'R': bc[1] = bc[4] = 1; break;
    case 'r': bc[1] = bc[4] = 1; break;
    case '1': bc[1] = bc[4] = 1; break;
    default : bc[1] = bc[4] = 0;
    }
    switch (*CONDLIM3) {
    case 'O': bc[2] = bc[5] = 2; break;
    case 'o': bc[2] = bc[5] = 2; break;
    case '2': bc[2] = bc[5] = 2; break;
    case 'R': bc[2] = bc[5] = 1; break;
    case 'r': bc[2] = bc[5] = 1; break;
    case '1': bc[2] = bc[5] = 1; break;
    default : bc[2] = bc[5] = 0;
    }
    sprintf (GRIDFILE, "ref_%s", ShortParName);
    if (!CPU_Rank || AllCPUs) {	/* TGF : temporary grid file */
      TGF = prs_openi (GRIDFILE);
      switch (NDIM) {
      case 1:
	fprintf (TGF, "%.15g\t%.15g\t%d %ld %ld\n",\
		 RANGE1LOW, RANGE1HIGH,\
		 0, bc[0], bc[3]);
	break;
      case 2:
	fprintf (TGF, "%.15g\t%.15g\t%.15g\t%.15g\t%d %ld %ld %ld %ld\n",\
		 RANGE1LOW, RANGE2LOW, RANGE1HIGH, RANGE2HIGH,\
		 0, bc[0], bc[1], bc[3], bc[4]);
	break;
      case 3:
	fprintf (TGF, "%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%d %ld %ld %ld %ld %ld %ld\n",\
		 RANGE1LOW, RANGE2LOW, RANGE3LOW, RANGE1HIGH, RANGE2HIGH, RANGE3HIGH,\
		 0, bc[0], bc[1], bc[2], bc[3], bc[4], bc[5]);
	break;
      }
      fclose (TGF);
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);	/* All processes will have to read the above written file */
  if ((NDIM < 1) || (NDIM > 3)) prs_error ("NDIM should be 1, 2 or 3.");
  if (NDIM < 3) INHIBREFDIM3 = TRUE;
  if (NDIM < 2) INHIBREFDIM2 = TRUE;
  Refine[0] = (INHIBREFDIM1 == TRUE ? FALSE : TRUE);
  Refine[1] = (INHIBREFDIM2 == TRUE ? FALSE : TRUE);
  Refine[2] = (INHIBREFDIM3 == TRUE ? FALSE : TRUE);
  for (i = 0; i < 3; i++)
    Nghost[i] = (i < NDIM ? NGH : 0);
  if (strncasecmp(ADIABATIC, "Y", 1) == 0) {
    Isothermal = NO;
    pInfo ("System is adiabatic\n");
    TOGPO = 2.0/(GAMMA+1.0);
    TOGMO = 2.0/(GAMMA-1.0);
    GMOGPO = (GAMMA-1.0)/(GAMMA+1.0);
    OOG = 1.0/GAMMA;
    OOGMO = 1.0/(GAMMA-1.0);
  }else{
    Isothermal = YES;
    pInfo ("System is (locally) isothermal\n");
  }

  if(Isothermal)
    __Compute_Fluxes = &Compute_Fluxes_Iso;
  else
    __Compute_Fluxes = &Compute_Fluxes_Adi;

  __ExecCommUpVar = &ExecCommUpVarLIL;
  if ((GHOSTFILLINGORDER > 1) || (GHOSTFILLINGORDER < 0))
    pWarning ("GHOSTFILLINGORDER can only be 1 (default) or 0 (direct injection). Set to 1.\n");
  if (GHOSTFILLINGORDER == 0)
    __ExecCommUpVar = &ExecCommUpVar;

  set = NO;
  if (strncasecmp(RIEMANNSOLVER, "2S", 2) == 0) {
    if(Isothermal){
      pInfo ("Two shock isothermal Riemann solver\n");
      __Riemann_Solver = &GetStar_TWOSHOCKS;
    }else{
      pInfo ("Riemann solver cannot be selected in adiabatic version.\n");
    }	
    set = YES;
  }
  if (strncasecmp(RIEMANNSOLVER, "2R", 2) == 0) {
    if(Isothermal){
      pInfo ("Two rarefaction isothermal Riemann solver\n");
      __Riemann_Solver = &GetStar_TWORAREFACTIONS;
    }else{
      pInfo ("Riemann solver cannot be selected in adiabatic version.\n");
    }	
    set = YES;
  }
  if (strncasecmp(RIEMANNSOLVER, "IT", 2) == 0) {
    if(Isothermal){
      pInfo ("Iterative isothermal Riemann solver\n");
      __Riemann_Solver = &GetStar_ITERATIVE;
    }else{
      pInfo ("Riemann solver cannot be selected in adiabatic version.\n");
    }
    set = YES;
  }
  if(!set) prs_error ("Invalid Riemann Solver code.");

  if (strncasecmp(METHOD, "PLM", 3) == 0) {
    if(Isothermal){
      pInfo ("Piecewise linear method\n");
      __Prepare_Riemann_States = &plm;
    }else{
      pInfo ("Piecewise linear method adiabatic\n");
      __Prepare_Riemann_States = &plm_adiab;
    }
    mPLM = YES;
  } else {
    if (strncasecmp(METHOD, "GFO", 3) == 0) {
      if(Isothermal){
	pInfo ("Godunov First Order method\n");
	__Prepare_Riemann_States = &gfo;
      }else{
	pInfo ("Godunov First Order method adiabatic\n");
	__Prepare_Riemann_States = &gfo_adiab;
      }
      mGFO = YES;
    } else {
      if (strncasecmp(METHOD, "MUS", 3) == 0) {
	if(Isothermal){
	  pInfo ("MUSCL-Hancock method\n");
	  __Prepare_Riemann_States = &muscl;
	}else{
	  pInfo ("MUSCL-Hancock method for adiabatic setup\n");
	  __Prepare_Riemann_States = &muscl_adiab;
	}
	mMUSCL = YES;
      } else {
	prs_error ("Invalid Godunov method code.");
      }
    }
  }
  Periodic[0] = DIM1PERIODIC;
  Periodic[1] = DIM2PERIODIC;
  Periodic[2] = DIM3PERIODIC;
  for (i = 0; i < NDIM; i++)
    if (Periodic[i]) 
      pInfo ("System is periodic in %s\n", CoordNames[CoordType*3+InvCoordNb[i]]);
  Ncell0[0] = SIZE1;
  Ncell0[1] = SIZE2;
  Ncell0[2] = SIZE3;
  if (((SIZE2 > 1) || (SIZE3 > 1)) && (NDIM == 1)) {
    pWarning ("The run is declared as 1D (NDIM=1)");
    pWarning ("But either SIZE2 or SIZE3 is not set to 1");
  }
  if ((SIZE3 > 1) && (NDIM < 3)) {
    pWarning ("The run is not 3D (NDIM < 3)");
    pWarning ("and yet SIZE3 is not set to 1");
  }
  corner_min0[0] = RANGE1LOW;
  corner_min0[1] = RANGE2LOW;
  corner_min0[2] = RANGE3LOW;
  corner_max0[0] = RANGE1HIGH;
  corner_max0[1] = RANGE2HIGH;
  corner_max0[2] = RANGE3HIGH;
  GridFriction[0] = GRIDFRICTION1;
  GridFriction[1] = GRIDFRICTION2;
  GridFriction[2] = GRIDFRICTION3;
  if (KEPLERIAN && ((_RAD_ >= NDIM) || (_AZIM_ >= NDIM))) {
    pWarning ("KEPLERIAN flag is meaningless in that case since either the radius or the azimuth is not expressed. I set it to false");
    KEPLERIAN= FALSE;
  }
  if (!(*POTENTIALCODE == 0)) 
    EXTERNALPOTENTIAL = YES;
  if (EXTERNALPOTENTIAL)
    pInfo ("An external potential is applied.\n");
  else
    pInfo ("No external potential is applied.\n");      
  if (EXTERNALPOTENTIAL) {
    set = NO;
    DatePotentialConstant = atof (FIXEDPOTENTIAL);
    *FIXEDPOTENTIAL = (char)toupper ((int)*FIXEDPOTENTIAL);
    if (*FIXEDPOTENTIAL  == 'N') { /* like 'NO' or 'NEVER' */
      DatePotentialConstant = (real)NTOT*(real)DT*1e20;
      pInfo ("Potential is never assumed constant\n");
      set = YES;
    }
    if ((*FIXEDPOTENTIAL  == 'A') || (*FIXEDPOTENTIAL  == 'Y')) {/* like 'ALWAYS' or 'YES' */
      DatePotentialConstant = -(real)NTOT*(real)DT*1e20;
      pInfo ("Potential is assumed to be always constant\n");
      set = YES;
    }
    if (!set) pInfo ("Potential is assumed to be constant after date : %g\n", \
		     DatePotentialConstant);
  }
  /* Subcycling considerations */
  if (*SUBCYCLING == 'A') {
    SubCycling = AUTO;
  } else if (*SUBCYCLING == 'N') {
    SubCycling = NONE;
  } else if (*SUBCYCLING == 'F') {
    SubCycling = FULL;
  } else {
    prs_error ("SUBCYCLING string parameter should be one of the following : AUTO, NONE, FULL");
  }
  /*******************************************/
  /*******************************************/
  /* Hydrostatic equilibrium enforcement below */
  /*******************************************/
  /*******************************************/
  NDimHydroStat = 0;
  for (i = 0; i < 3; i++) {
    SqueezeDim[i] = (*(SYMDIM_EQ+i) == '1' ? TRUE: FALSE);
    CorrHydroStat[i] = (*(CORRDIM_EQ+i) == '1' ? TRUE: FALSE);
    if (CorrHydroStat[i]) {
      if (SqueezeDim[i]) {
	sprintf (message,"System is invariant in %s but you want to enforce hydrostatic equilibrium along this dimension.", \
		 CoordNames[CoordType*3+InvCoordNb[i]]);
	prs_error (message);
      }
      dimcorr[NDimHydroStat] = i;
      NDimHydroStat++;
      pInfo ("Hydrostatic equilibrium enforced in %s\n",\
	     CoordNames[CoordType*3+InvCoordNb[i]]);
      HydroStaticEnforce = YES;
    }
    if (SqueezeDim[i]) pInfo ("Hydrostatic system is invariant in %s\n", CoordNames[CoordType*3+InvCoordNb[i]]);
  }
  if (HydroStaticEnforce && (*POTENTIALHYDROSTAT == 0)) {
    pInfo ("Hydrostatic equilibrium is enforced but no potential law is defined by POTENTIALHYDROSTAT\n");
    pInfo ("I assume that it is the same law as defined by POTENTIALCODE\n");
    strcpy (POTENTIALHYDROSTAT, POTENTIALCODE);
  }
  if (HydroStaticEnforce && (*INITHYDROSTAT == 0)) {
    pInfo ("Hydrostatic equilibrium is enforced but no initial conditions are defined by INITHYDROSTAT\n");
    pInfo ("I assume that they are the same as defined by INITCODE\n");
    strcpy (INITHYDROSTAT, INITCODE);
  }
  /*******************************************/
  if (KEPLERIAN && !MERIDIANCORRECTION) {
    pInfo ("You have a Keplerian disk, but you have deactivated\n");
    pInfo ("the meridian correction that yields a wakeless\n");
    pInfo ("system of nested meshes\n");
  }
  /*******************************************/
  /*******************************************/
  if ((!MERIDIANCORRECTION) && (!KEPLERIAN)) {
    pInfo ("You have deactivated the meridian fluxes corrections but this would have\n");
    pInfo ("an effect only in a Keplerian disk (KEPLERIAN boolean set to True)\n");
  }
  if ((KEPLERIAN)  && (!TorqueMon)) {
    pInfo ("\nYou have conditions specific of a Keplerian disk\n");
    pInfo ("but you have not activated the torque monitoring (-t)\n");
    pInfo ("I assume you wanted to activate it.\n\n");
    TorqueMon = TRUE;
  }
  MultiFluid ();       	/* Check how many fluids */
  /* The call to the above function must be performed *after* the
     examination of INITCODE/INITHYDROSTAT variables*/
  *NODAMPING = (char)toupper ((int)*NODAMPING);
  if( (*NODAMPING == 'Y')  || (*NODAMPING == 'T') || (*NODAMPING == '1'))
    NoStockholm = YES;
  if (NoStockholm == YES) {
    if (KEPLERIAN) {
      pInfo ("\nAlthough you have conditions specific of a Keplerian disk,\n");
      pInfo ("you have deactivated the damping boundary conditions in radius\n");
      pInfo ("and colatitude.\n\n");
    } else {
      pInfo ("\nYou have set the '-h' flag to deactivate damping boundary conditions.\n");
      pInfo ("This is however irrelevant since the KEPLERIAN flag is set to false\n");
      pInfo ("so no damping boundary conditions would be applied anyhow.\n\n");
    }
  } else {
    if (KEPLERIAN) {
      pInfo ("\nKEPLERIAN flag is set.\n");
      pInfo ("Damping boundary conditions will be applied to the disk\n\n");
    }
  }
  if ((strncmp (RADIATIVE, "YES", 3) == 0)  || (strncmp (RADIATIVE, "Yes", 3) == 0) || (strncmp (RADIATIVE, "True", 3) == 0) || (strncmp(RADIATIVE, "TRUE", 3) == 0)) Radiative = YES;
  if ((strncmp (STELLAR, "YES", 3) == 0)  || (strncmp (STELLAR, "Yes", 3) == 0) || (strncmp (STELLAR, "True", 3) == 0) || (strncmp(STELLAR, "TRUE", 3) == 0)) Stellar = YES;
  if ((strncmp (HALFDISK, "YES", 3) == 0)  || (strncmp (HALFDISK, "Yes", 3) == 0) || (strncmp (HALFDISK, "True", 3) == 0) || (strncmp(HALFDISK, "TRUE", 3) == 0)) HalfDisk = YES;
  if ((!Restart) && SuperImpose)
    prs_error ("You need to do a restart to activate superimposition.");
  return;}

