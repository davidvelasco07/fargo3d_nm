//MPI global variables

int NCPUS_per_domain;
int NFluids_per_rank;
int NDomains;

int CPU_World_Rank;
int CPU_World_Number;

int CPU_Rank;    //For each fluid (rows)
int CPU_Number;  //For each fluid (rows)

int CPU_RowRank;   //For each fluid (rows)
int CPU_RowNumber; //For each fluid (rows)

int CPU_ColumnRank;    //For each domain (columns)
int CPU_ColumnNumber;  //For each domain (columns)

MPI_Comm DomainComm; //Intra communicator for each fluid
MPI_Comm FluidsComm; //Intra communicator for each domain

MPI_Request RequestTotalDensity;
MPI_Request RequestGasDensity;
MPI_Request RequestGasEnergy;
MPI_Request request_c;

int FluidColor;   //Rows
int DomainColor;  //Columns

int Current_Level = 0;
int Current_Fluid = 0;
tGrid_CPU *Current_Grid;
tGrid_CPU *Grid_item;
int NbRestart = 0;
boolean CPU_Master = YES;
boolean CPU_Global_Master = YES;
int nbdump;

//Global variables
boolean FirstPassXAdvection = NO;
boolean Resistivity_Profiles_Filled = NO;
boolean VxIsResidual = NO;
boolean LogGrid = NO;
boolean GuidingCenter = NO;
boolean Corotating = NO;
boolean Restart = NO;
boolean Restart_Full = NO;
boolean Stockholm = NO;
boolean Merge = YES;
boolean Merge_All = NO;
boolean MonitorIntegral;
boolean TimeInfo = NO;
boolean EverythingOnCPU = NO;
boolean ForwardOneStep = NO;
boolean Vtk2dat = NO;
boolean Dat2vtk = NO;
boolean PostRestart = NO;
boolean OnlyInit = NO;
boolean EarlyOutputRename = NO;
boolean RedefineOptions = NO;
boolean DeviceFileSpecified = NO;
boolean StretchOldOutput = NO;
boolean ThereArePlanets = NO;
boolean ThereIsACentralBinary = NO;
boolean AddSubPatch = NO;
real    PhysicalTimeInitial;
real    PhysicalTime = 0;
real    XAxisRotationAngle = 0.0;
char    NewOutputdir[1024];
char    VersionString[1024];
char    StickyOptions[1024];
char    BoundaryFile[4096];
char    CommandLine[1024];
char    FirstCommand[1024];
char    ArchFile[1024];
char    ParameterFile[1024];
char    DefaultOut[1024];
char    DeviceFile[1024];
char    CurrentWorkingDirectory[1024];
char    *InputFile;
char    *PlanetaryFile;
char    *HostsList;
int     TimeStep = 0;
int     FullArrayComms = 0;
int     ContourComms = 0;
int     DeviceManualSelection = -1;
int     ArrayNb = 0;
int     StretchNumber = 0;
int     BinaryStar1 = 0;
int     BinaryStar2 = 0;

real Domega;

real InnerBorder;
real OuterBorder;

TimeProcess t_speedup_cpu;
TimeProcess t_speedup_gpu;
int t_speedup_count;
real time_speedup_cpu;
real time_speedup_gpu;

PlanetarySystem *Sys;
Point DiskOnPrimaryAcceleration;
Point IndirectTerm;

real StepTime;

real localforce[12];
real globalforce[12];

//FIELDS VARIABLES

Grid Gridd;

real *Fluxes[3][2];
Field *Flux = NULL;

Field *Vx=NULL;
Field *Vy=NULL;
Field *Vz=NULL;

Field *Vx_temp=NULL;
Field *Vy_temp=NULL;
Field *Vz_temp=NULL;

Field *Vx_half=NULL;
Field *Vy_half=NULL;
Field *Vz_half=NULL;

Field *Slope=NULL;

Field *Mpx=NULL;
Field *Mpy=NULL;
Field *Mpz=NULL;
Field *Mmx=NULL;
Field *Mmy=NULL;
Field *Mmz=NULL;

Field *Pot=NULL;

Field *DivRho=NULL;
Field *DensStar=NULL;
Field *Qs=NULL;

Field *Density=NULL;
Field *Energy=NULL;
Field *Pressure=NULL;

Field *Total_Density=NULL;

Field *QL=NULL;
Field *QR=NULL;
Field *QLE=NULL;
Field *QRE=NULL;
Field *LapPPA=NULL;

Field *Sdiffyczc=NULL;
Field *Sdiffyfzc=NULL;
Field *Sdiffyczf=NULL;
Field *Sdiffyfzf=NULL;

Field *ThermalDiff=NULL;

// Below: fields specific to FARGO algorithms
Field2D *VxMed=NULL;
Field2D *Vxhy=NULL;
Field2D *Vxhyr=NULL;
Field2D *Vxhz=NULL;
Field2D *Vxhzr=NULL;
Field2D *Reduction2D=NULL;

FieldInt2D *Nxhy=NULL;
FieldInt2D *Nxhz=NULL;
FieldInt2D *Nshift=NULL;

//MHD FIELDS
//#ifdef MHD
Field *Bx=NULL;
Field *By=NULL;
Field *Bz=NULL;

Field *B1_star=NULL;
Field *B2_star=NULL;

Field *V1_star=NULL;
Field *V2_star=NULL;

Field *Slope_b1=NULL;
Field *Slope_v1=NULL;
Field *Slope_b2=NULL;
Field *Slope_v2=NULL;

Field *Emfx=NULL; 
Field *Emfy=NULL;
Field *Emfz=NULL;

Field *EmfxH=NULL;
Field *EmfyH=NULL;
Field *EmfzH=NULL;
Field *BxH=NULL;
Field *ByH=NULL;
Field *BzH=NULL;
Field *Jx=NULL;
Field *Jy=NULL;
Field *Jz=NULL;
Field *EtaOhm=NULL;
Field *EtaHall=NULL;
Field *EtaAD=NULL;

Field *Divergence=NULL;
//#endif

Field2D *Density0=NULL;
Field2D *Vx0=NULL;
Field2D *Vy0=NULL;
Field2D *Vz0=NULL;
Field2D *Energy0=NULL;

//Communications variables

int Ncpu_x; // Numbers of cpus in x-axis;
int Ncpu_y; // Idem in y-axis;

Buffer Bfd;   //|
Buffer Bfu;   //|
Buffer Bfl;   //|
Buffer Bfr;   //|----> interface buffers for density field; (borders)
Buffer Bfcdl; //|                view comm.c
Buffer Bfcdr; //|
Buffer Bfcul; //|
Buffer Bfcur; //|

//Useful numbers

//CPU GLOBAL LIGHT ARRAYS

real Dx;
real *Xmin=NULL;
real *Ymin=NULL;
real *Zmin=NULL;
real *Xmed=NULL;
real *Ymed=NULL;
real *Zmed=NULL;
real *InvDiffXmed=NULL;
real *InvDiffYmed=NULL;
real *InvDiffZmed=NULL;
real *Sxj=NULL;
real *Sxk=NULL;
real *Syj=NULL;
real *Syk=NULL;
real *Szj=NULL;
real *Szk=NULL;
real *InvVj=NULL;
real shift_buffer[MAX1D];

//GPU GLOBAL LIGHT ARRAYS
real *Alpha=NULL;
real *Alpha_d=NULL;
real *Dx_d=NULL;
real *Xmin_d=NULL;
real *Ymin_d=NULL;
real *Zmin_d=NULL;
real *Sxj_d=NULL;
real *Sxk_d=NULL;
real *Syj_d=NULL;
real *Syk_d=NULL;
real *Szj_d=NULL;
real *Szk_d=NULL;
real *InvVj_d=NULL;
real shift_buffer_d[MAX1D];

//Grid variables

int Nx;
int Ny;
int Nz;
int J;
int K;
int Y0;
int Z0;
int Stride_cpu;
int Stride_gpu;
int Pitch_cpu;
int Pitch_gpu;
int Pitch_Int_gpu;
int Pitch2D;
int Stride;
int ix;
int iy;
int ycells;
int zcells;
int y0cell;
int z0cell;
int Maxsize_cpu = 1;
int Maxsize2D_cpu = 1;
int Maxsize_gpu = 1;
int Maxsize2D_gpu = 1;
int MaxsizeInt_cpu = 1;
int MaxsizeInt_gpu =1;
int LocalBC[3][2];

//Entry point to JUPITER architecture
tGrid_CPU *Current_Jupiter_Patch = NULL;

//For checknan
Field *ListOfGrids = NULL;

//Boundary variables

int Bounl;
int Bounr;
int Bounu;
int Bound;

//psys variables.
real Xplanet;
real Yplanet;
real Zplanet;
real VXplanet;
real VYplanet;
real VZplanet;
real MplanetVirtual;

MPI_Status fargostat;

real OMEGAFRAME0;

int Fscan;

long VtkPosition = 0; 

//Multifluid variables
int Timestepcount = 0;
int Fluidtype;
int FluidIndex;
real *Coeffval;
real Min[NFLUIDS];      //Comment: NFLUIDS is the upper bound for the size of the array.
Fluid *Fluids[NFLUIDS]; //Comment: NFLUIDS is the upper bound for the size of the array.
real M_acc[NFLUIDS][100];
real M_dot[NFLUIDS][100];

//Nested Meshes variables
char CoordNames[9][80] = {"X", "Y", "Z", "Radius", "Azimuth", "Z", "Radius", "Azimuth", "Co-latitude"};
char SCoordNames[9][80]= {"X", "Y", "Z", "Rad", "Phi", "Z", "Rad", "Phi", "Theta"};
char OutputDir[MAXLINELENGTH];
char SubPatchInfo[MAXLINELENGTH];
int NDIM;
int InvCoordNb[3];
int CoordNb[3];
int CoordType;
long MaxLevel = 100;
real corner_min0[3];
real corner_max0[3];
long ncorner_min0[3];
long ncorner_max0[3];
long jstride[3];
long LevMax = 0;
boolean Refine[3];
boolean Periodic[3];
long Ncell0[3];
long Nghost[3];
boolean AllowFlushLog = YES;
tGrid *GridList = NULL;
tGrid_CPU *Grid_CPU_list = NULL;
int _X_;
int _Y_;
int _Z_;
int _RAD_;
int _AZIM_;
int _COLAT_;
int _Radial_;
int _Azimuthal_;
int _Vertical_;
int _Colatitude_;
long cpugrid_number;
jCommunicator *ComListGhost = NULL;
jCommunicator *ComListFlux = NULL;
jCommunicator *ComListMean = NULL;
CommHash **CommHashSrc;  
int NbFluids = 1;
long TimeStepRatio[100];
long FluxIndex = 0;
real MaxLowLevDT[100];
boolean LevelHasChangedSinceCFL[100];
real MaxLowLevDT_loc[100];
long TimeStepRatio[100];
long BaseStepRatio[100];
long LevelCost[100];
long SubCycling;
real LevelDate[100];
long centered[80][3];
int staggered[80];
int stdim[80];
real *source[80];
real *dest[80];
int Ngrids=0;

//Pointers to functions
//WARNING!!! FUNCTIONS' ARGUMENTS MUST NOT CONTAIN BLANK SPACES
void (*ComputePressureFieldIso)();
void (*ComputePressureFieldAd)();
void (*ComputePressureFieldPoly)();
void (*SubStep1_x)(real);
void (*SubStep1_y)(real);
void (*SubStep1_z)(real);
void (*SubStep2_a)(real);
void (*SubStep2_b)(real);
void (*SubStep3)(real);
void (*SubStep4_a)(real);
void (*SubStep4_b)(real);
void (*Compute_ThermalDiffusion)();
void (*DivideByRho)(Field*,Field*);
void (*VanLeerX_a)(Field*);
void (*VanLeerX_b)(real,Field*,Field*,Field*);
void (*VanLeerY_a)(Field*);
void (*VanLeerY_b)(real,Field*,Field*);
void (*VanLeerZ_a)(Field*);
void (*VanLeerZ_b)(real,Field*,Field*);
void (*momenta_x)();
void (*momenta_y)();
void (*momenta_z)();
void (*reduction_SUM)(Field*,int,int,int,int); 
void (*reduction_MIN)(Field*,int,int,int,int); 
void (*Update_a_X)(real,Field*,Field*,Field*);
void (*Update_dens_a_X)(real,Field*,Field*,Field*);
void (*Update_b_X)(Field*,Field*);
void (*Update_a_Y)(real,Field*,Field*);
void (*Update_dens_a_Y)(real,Field*,Field*);
void (*Update_b_Y)(Field*,Field*);
void (*Update_a_Z)(real,Field*,Field*);
void (*Update_dens_a_Z)(real,Field*,Field*);
void (*Update_b_Z)(Field*,Field*);
void (*NewVelocity_x)();
void (*NewVelocity_y)();
void (*NewVelocity_z)();
void (*AdvectSHIFT)(Field*,FieldInt2D*);
void (*ComputeResidual)(real);
void (*ChangeFrame)(int,Field*,Field2D*);
void (*Potential)();
void (*CorrectVtheta)(real);
void (*cfl)(void);
void (*_ComputeForce)(real,real,real,real,real,char*);
void (*copy_velocities)(int);
void (*VanLeerX_PPA_a)(Field*);
void (*VanLeerX_PPA_b)(Field*);
void (*VanLeerX_PPA_steep)(Field*);
void (*VanLeerX_PPA_c)(Field*);
void (*VanLeerX_PPA_d)(real,Field*,Field*,Field*);
void (*VanLeerX_PPA_d_2d)(real,Field*,Field*,Field2D*);
void (*mon_dens)();
void (*mon_momx)();
void (*mon_momy)();
void (*mon_momz)();
void (*mon_torq)();
void (*mon_reynolds)();
void (*mon_maxwell)();
void (*mon_bxflux)();
void (*comm)();
void (*Reset_field)(Field*);
void (*ComputeTotalDensity)();
void (*copy_field)(Field*,Field*);
//DUST DIFFUSION
void (*DustDiffusion_Core)(real);
void (*DustDiffusion_Coefficients)();
//MHD..........................................
void (*ComputeSlopes)(int,int,int,Field*,Field*);
void (*_ComputeStar)(real,int,int,int,int,int,int,int,int,int,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*);
void (*_ComputeEmf)(real,int,int,int,int,int,int,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*,Field*);
void (*_UpdateMagneticField)(real,int,int,int,int,int,int,int,int,int,Field*,Field*,Field*);
void (*_LorentzForce)(real,int,int,int,int,int,int,int,int,int,int,int,Field*,Field*,Field*,Field*,Field*);
void (*_Resist)(int,int,int,int,int,int,int,int,int,Field*,Field*,Field*,Field2D*);
void (*EMF_Upstream_Integrate)(real);

void (*ComputeJx)();
void (*ComputeJy)();
void (*ComputeJz)();

void (*_OhmicDiffusion_emf)(int,int,int,int,int,int,int,int,int,Field*,Field*,Field*);
void (*OhmicDiffusion_coeff)();
void (*HallEffect_emfx)();
void (*HallEffect_emfy)();
void (*HallEffect_emfz)();
void (*_HallEffect_UpdateB)(real,int,int,int,int,int,int,int,int,int,int,int,Field*,Field*,Field*);
void (*HallEffect_UpdateEmfs)();
void (*HallEffect_coeff)();
void (*AmbipolarDiffusion_emfx)();
void (*AmbipolarDiffusion_emfy)();
void (*AmbipolarDiffusion_emfz)();
void (*AmbipolarDiffusion_coeff)();
//.............................................

void (*StockholmBoundary)(real);
void (*visctensor_cart)();
void (*addviscosity_cart)(real);
void (*visctensor_cyl)();
void (*addviscosity_cyl)(real);
void (*visctensor_sph)();
void (*addviscosity_sph)(real);
void (*Fill_GhostsX)();
void (*CheckMuteY)();
void (*CheckMuteZ)();

void (*SetupHook1)();

// Collisions
void (*_collisions)(real,int,int,int,int);
void (*ComputeDragCoeff)();
void (*ComputeCBcollisions_c)(real);
void (*_ComputeCBcollisions_cv)(real,int,int,int,Field*,Field*);
void (*_UpdateVelcollisions)(real,int,int,int,Field*,Field*);

//DragForce
void (*DragForce_Coeff)();
void (*DragForce_SumC)(real);
void (*_DragForce_SumCV)(real,int,int,int,Field*,Field*);
void (*_DragForce_UpdateVel)(real,int,int,int,Field*,Field*);


void (*Floor)();

void (*__WriteField)();
void (*__Restart)(Field*,int);

//Comment: NFLUIDS is the upper bound for the size of the array.
void (*boundary_xmin[NFLUIDS])();
void (*boundary_xmax[NFLUIDS])();
void (*boundary_ymin[NFLUIDS])();
void (*boundary_ymax[NFLUIDS])();
void (*boundary_zmin[NFLUIDS])();
void (*boundary_zmax[NFLUIDS])();

void (*Accretion)(real,real,real);
void (*PlanetHeating)(real,real,real,real);

int Grid_id;