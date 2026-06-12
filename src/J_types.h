//#define MAXLINELENGTH 512L
//#define MAXNAMELENGTH 128L
//
//
//typedef double real;
//typedef long   boolean;
//
//struct param {
//  char name[256];
//  long type;
//  char *variable;
//  long read;
//  long necessary;
//};
//
struct pair {
  real min;
  real max;
};

//
//typedef struct param Param;
//
//struct hashparam {
//  char name[MAXLINELENGTH];
//  char stringvalue[MAXLINELENGTH];
//  real floatvalue;
//  long intvalue;
//  long boolvalue;
//};
//
//typedef struct hashparam HashParam;
//
//struct radixstruct {
//	real	constant;
//	real 	radius_slope;
//	real	sintheta_slope;
//};
//typedef struct radixstruct RadixStruct;
//
enum {_density_, _energy_, _vrad_, _vazimuth_, _vcolatitude_, _other_, _tau_, _erad_, _qp_};
//
enum {_Density_, _Energy_, _Vx_, _Vy_, _Vz_, _Potential_, _Erad_, _QP_, _EradDeriv_,_Label_};
///* Note: if you add additional fields in the above enum, check straight away that 
//the Ptr[] array of FluidPatch is large enough (see gridtypes.h) */
