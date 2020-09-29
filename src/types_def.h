typedef char boolean;
#ifdef FLOAT
typedef float real;
#else
typedef double real;
#endif
typedef struct param Param;
typedef struct planetary_system PlanetarySystem;
typedef struct field Field;
typedef struct field2D Field2D;
typedef struct fieldint2D FieldInt2D;
typedef struct grid Grid;
typedef struct buffer Buffer;
typedef struct force Force;
typedef struct point Point;
typedef struct timeprocess TimeProcess;
typedef struct zmeanprop ZmeanProp;
typedef struct hashparam HashParam;
typedef struct state_vector StateVector;
typedef struct orbital_elements OrbitalElements;

typedef struct fluid Fluid;

typedef struct gridfileinfo GridFileInfo;
typedef struct fluidpatch FluidPatch;

typedef struct tgrid tGrid;
typedef struct tgrid_cpu tGrid_CPU;
typedef struct pair Pair;
typedef struct destblocks DestBlocks;
typedef struct scalarfield ScalarField;
typedef struct vectorfield VectorField;
typedef struct jcommunicator jCommunicator;
typedef struct commhash CommHash;

struct pair {
  real min;
  real max;
};

enum {_density_, _energy_, _vrad_, _vazimuth_, _vcolatitude_, _other_, _tau_, _erad_, _qp_};
//
enum {_Density_, _Energy_, _Vx_, _Vy_, _Vz_, _Vx_temp_, _Vy_temp_, _Vz_temp_,_Label_};
///* Note: if you add additional fields in the above enum, check straight away that 
//the Ptr[] array of FluidPatch is large enough (see gridtypes.h) */