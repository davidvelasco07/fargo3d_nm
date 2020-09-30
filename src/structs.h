struct param{
  char name[80];
  char *variable;
  int type;
  int need;
  int read;
};


struct hashparam {
  char name[MAXLINELENGTH];
  char stringvalue[MAXLINELENGTH];
  real floatvalue;
  long intvalue;
  long boolvalue;
};



/** Contains all the information about a planetary system at a given
    instant in time. */
struct planetary_system {
  int nb;			/**< Number of planets */
  real *mass;			/**< Masses of the planets */
  real *mass_cpu;			/**< Masses of the planets */
  real *mass_gpu;			/**< Masses of the planets */
  real *x;			/**< x-coordinate of the planets */
  real *x_cpu;			/**< x-coordinate of the planets */
  real *x_gpu;			/**< x-coordinate of the planets */
  real *y;			/**< y-coordinate of the planets */
  real *y_cpu;			/**< y-coordinate of the planets */
  real *y_gpu;			/**< y-coordinate of the planets */
  real *z;			/**< z-coordinate of the planets */
  real *z_cpu;			/**< z-coordinate of the planets */
  real *z_gpu;			/**< z-coordinate of the planets */
  real *vx;			/**< x-coordinate of the planets'velocities */
  real *vy;		        /**< y-coordinate of the planets'velocities */
  real *vz;		        /**< z-coordinate of the planets'velocities */
  real *acc;			/**< The planets' accretion times^-1 */
  char **name;  		/**< The planets' names */
  boolean *FeelDisk;		/**< For each planet tells if it feels the disk (ie migrates) */
  boolean *FeelOthers;		/**< For each planet tells if it feels
				   the other planets' gravity */
};

struct grid {   // Store all relevant information for grid
  int nx; //Nsec --> j
  int ny; //Nrad --> i
  int nz; //Ncol --> k
  int J;  //|--> Global index in X
  int K;  //| --> Global index in Y
  int stride; // init in CreateField();
  int bc_up;
  int bc_down;
  int bc_left;
  int bc_right;
  int NJ; // --> Size of array of processors, in Y (same for everyone)
  int NK; // --> Size of array of processors, in Z (same for everyone)
};

struct zmeanprop {
  real PreviousDate;
  real ResetDate;
  boolean NeverReset;
};

struct field { //Multiple fields on code (density, vx, vy,...)
  char *name;
  real *field_cpu;
  real *backup; // used for creation of check points in debugging CPU vs GPU
  real *secondary_backup; // same thing
  real *field_gpu;
  int x_cpus;
  int y_cpus;
  int type;
  int level;
  boolean *fresh_cpu;
  boolean *fresh_gpu;
  boolean fresh_inside_contour_cpu[4];
  boolean fresh_inside_contour_gpu[4];
  boolean fresh_outside_contour_cpu[4];
  boolean fresh_outside_contour_gpu[4];
  struct field *next; // Linkedlist
  struct field **owner; //used for aliases
#ifdef GPU
  struct cudaPitchedPtr gpu_pp;
  struct cudaPitchedPtr cpu_pp;
#endif
  char file_origin[MAXLINELENGTH];
  int line_origin;
  real *x;
  real *y;
  real *z;
  double *zmean;
  struct zmeanprop zp;
  tGrid_CPU *desc;
};

struct fluid {
  char *name;
  int Fluidtype;
  struct field2D *VxMed;
  struct field *Density;
  struct field *Energy;
  struct field *Vx;
  struct field *Vx_temp;
  struct field *Vy;
  struct field *Vy_temp;
  struct field *Vz;
  struct field *Vz_temp;
  struct field *Vx_half;
  struct field *Vy_half;
  struct field *Vz_half;
#ifdef STOCKHOLM
  struct field2D *Density0;
  struct field2D *Energy0;
  struct field2D *Vx0;
  struct field2D *Vy0;
  struct field2D *Vz0;
#endif
};

struct field2D { //Multiple 2D fields on code (azimuthal averages, etc.)
  char *name;
  real *field_cpu;
  real *backup; // used for creation of check points in debugging CPU vs GPU
  real *secondary_backup;
  real *field_gpu;
  int x_cpus;
  int y_cpus;
  size_t pitch;
  int kind;
  boolean *fresh_cpu;
  boolean *fresh_gpu;
  struct field2D *next; // Linkedlist
  long size;
  tGrid_CPU *desc;
};

struct fieldint2D { //Multiple 2D fields on code (azimuthal averages, etc.)
  char *name;
  int *field_cpu;
  int *backup;
  int *secondary_backup;
  int *field_gpu;
  int x_cpus;
  int y_cpus;
  boolean *fresh_cpu;
  boolean *fresh_gpu;
  struct fieldint2D *next; // Linkedlist
};

struct buffer {
  real *buffer;
  int index;
};

struct force {			
  real fx_inner;		/**< x-component of the force arising from the inner disk, without Hill sphere avoidance  */
  real fy_inner;		/**< y-component of the force arising from the inner disk, without Hill sphere avoidance  */
  real fz_inner;
  real fx_ex_inner;    /**< x-component of the force arising from the inner disk, with Hill sphere avoidance  */	
  real fy_ex_inner;    /**< y-component of the force arising from the inner disk, with Hill sphere avoidance  */        
  real fz_ex_inner;
  real fx_outer;         /**< x-component of the force arising from the outer disk, without Hill sphere avoidance */
  real fy_outer;	        /**< y-component of the force arising from the outer disk, without Hill sphere avoidance */
  real fz_outer;
  real fx_ex_outer;	/**< x-component of the force arising from the outer disk, with Hill sphere avoidance  */	
  real fy_ex_outer;	/**< x-component of the force arising from the outer disk, with Hill sphere avoidance  */        
  real fz_ex_outer;
};


struct point {
  real x;
  real y;
  real z;
};

/** This structure is used for monitoring CPU time usage. It is used
    only if -t is specified on the command line. */
struct timeprocess {
  char name[80];
  clock_t clicks;
};

struct orbital_elements {
  real a;
  real e;
  real i;
  real an;
  real per;
  real ta;
  real E; //Additional useful anomalies
  real M;
  real Perihelion_Phi;
};

struct state_vector {
  real x;
  real y;
  real z;
  real vx;
  real vy;
  real vz;
};

/** The mesh structure that contains the information independent on
    the number of processing elements */
struct gpucomm{
  int facedim;
  int faceside;
  int size;
  int sizeD;
  int pitchD;
  int strideD;
  int Pitch;
  int Stride;
  int xmin;
  int ymin;
  int zmin;
  int cxmin;
  int cymin;
  int czmin;
  int dx;
  int dy;
  int dz;
  int srcg;
};

struct gpucommDown{
  int facedim;
  int faceside;
  int size;
  int sizeD;
  int pitchD;
  int strideD;
  int Pitch;
  int Stride;
  int xmin;
  int ymin;
  int zmin;
  int xrefine;
  int yrefine;
  int zrefine;
  int srcxface;
  int srcyface;
  int srczface;
  int srcg;
};

struct gpucommFlux{
  int facedim;
  int faceside;
  int size;
  int sizeD;
  int pitchD;
  int strideD;
  int n1;
  int xmin;
  int ymin;
  int zmin;
  int xrefine;
  int yrefine;
  int zrefine;
  int sizeS;
  int srcg;
};

struct parms{
  int pitch;
  int stride;
  int cxmin;
  int cymin;
  int czmin;
  int dx;
  int dy;
  int dz;
  int sizex;
  int sizey;
  int sizez;
};

struct gridcomms{
  struct gpucomm *ParmsUp;
  struct gpucommDown *ParmsDown;
  struct gpucommFlux *ParmsFlux;
  struct commhash *CommListSame;
  struct commhash *CommListUp;
  struct commhash *CommListDown;
  struct commhash *CommListFlux;
  real **BuffersUp;
  real **BuffersDown;
  real **BuffersFlux;
  int nbCommsUp;
  int nbCommsDown;
  int nbCommsFlux;
};

struct tgrid {
  long level;			/**< Refinement level */
  long ncell[3];		/**< Number of zones in each dimension (excluding ghosts) */
  long gncell[3];		/**< Number of zones in each dimension (including ghosts) */
  long nsize[3];		/**< Absolute size (excluding ghosts) */
  long ncorner_min[3];		/**< Absolute position (min corner, excluding ghosts) */
  long ncorner_max[3];		/**< Absolute position (max corner, excluding ghosts) */
  long gnsize[3];		/**< Absolute size (including ghosts) */
  long gncorner_min[3];		/**< Absolute position (min corner, including ghosts) */
  long gncorner_max[3];		/**< Absolute position (max corner, including ghosts) */
  real size[3];			/**< Size in physical units */
  real corner_min[3];		/**< Position of min corner in physical units */
  real corner_max[3];		/**< Position of max corner in physical units */
  real *Edges[3];		/**< Interface position in each dimension */
  boolean present;		/**< True if and only if that grid is mapped */
				/**< by local CPU */
  long number;			/**< Logical number of grid */
  struct tgrid *next;		/**< Pointer to next grid */
  struct tgrid *prev;		/**< Pointer to previous grid */
  long BoundaryConditions[3][2];/**< List of BC's on faces */
  long minCPU;			/**< range of processes over which */
  long maxCPU;			/**< the grid is distributed */
  long dn[3];			/**< Cell size in absolute units */
  long linenumber;		/**< Line number in grid file */
  long monoCPU;			/**< Information used for a restart only */
  long Ncpus[3];  /** Number of subdomains per direction */
};

/** The information relative to a fluid patch that is accessible to a
given processing element.  The 'parent' (or 'Parent') information
refers to the tGrid of which the current FluidPatch is the local
subset. */

struct tgrid_cpu {
  long   parent;       		/**< Logical number of "parent" grid */
  long   number;		/**< Internal logical number */
  tGrid *Parent;		/**< Pointer to "parent */
  long   level;			/**< Refinement level */
  int    cpu;			/**< process number of corresponding CPU grid */
  long   ncell[3];		/**< Number of zones in each dimension (excluding ghosts)*/
  long   gncell[3];		/**< Number of zones in each dimension (including ghosts)*/
  long   stride[3];		/**< Data storage stride in each dimension */
  long   nsize[3];		/**< Absolute size */
  long dn[3];			/**< Cell size in absolute units */
  long iface[3][2];		/**< Interface type (CPU->brother, True BC, Mesh->other mesh) */
  long ncorner_min[3];		/**< Absolute position (min corner) */
  long ncorner_max[3];		/**< Absolute position (max corner) */
  long pcorner_min[3];		/**< Position in parent (in cells) - min */
  long pcorner_max[3];		/**< Position in parent (in cells) - max */
  long gnsize[3];		/**< Absolute size (including ghosts) */
  long gncorner_min[3];		/**< Absolute position (min corner, including ghosts) */
  long gncorner_max[3];		/**< Absolute position (max corner, including ghosts) */
  real size[3];			/**< Size in physical units */
  real corner_min[3];		/**< Position of min corner in physical units */
  real corner_max[3];		/**< Position of max corner in physical units */
  real *InvVolume;		/**< Inverse of volume of each cell */
  real *Center[3];		/**< Center of zones in each dimension */
  real *Metric[3][2];		/**< 1D Arrays for metric coefficients */
  real *InvMetric[3][2];	/**< 1D Arrays for metric coefficients */
  real *InterSurface[3];	/**< Surface of interfaces in each dimension */
  real *Edges[3];		/**< Interface position in each dimension */
  char *Hidden;			/**< A set of flags which says if the zone lies behind another */
  FluidPatch *fluid;		/**< Pointer to the first fluid patch (follow with ->next for multifluid) */
  struct fluid *Fluids[NFLUIDS];
  int color;                    /**< Variables required for the Stellar Irradiation with MPI */
  int colorz;                   /**< Variables required for the Stellar Irradiation with MPI */
  int key;                      /**< Variables required for the Stellar Irradiation with MPI */
  real *optical_depth;          /**< Total optical depth of the patch at stellar wavelengths along radial direction */ 
  struct tgrid_cpu *next;	/**< Pointer to next CPU Grid */
  struct tgrid_cpu *prev;	/**< Pointer to previous CPU Grid */
  int Pitch_gpu;
  int Stride_gpu;
  int Pitch_cpu;
  int Stride_cpu;
  int Pitch2D;
  int PitchInt;
  real *Xmed;
  real *Ymed;
  real *Zmed;
  real *InvDiffXmed;
  real *InvDiffYmed;
  real *InvDiffZmed;
  real *Sxj;
  real *Sxk;
  real *Syj;
  real *Syk;
  real *Szj;
  real *Szk;
  real *InvVj;
  int LumCells;
#ifdef GPU
  struct gridcomms src;
  struct gridcomms dst;
  struct parms gpuparms;
  real **source;
  real **dest;
  int *centered;
  int *centered_gpu;
  real *Xmin_d;
  real *Ymin_d;
  real *Zmin_d;
  real *Sxj_d;
  real *Sxk_d;
  real *Syj_d;
  real *Syk_d;
  real *Szj_d;
  real *Szk_d;
  real *InvVj_d;
  char *Hidden_d;
#endif
};

/** Container used to synchronize the ghost zones between different
grids.  A given communicator is used by two different CPU_Grids (and
two only), and it contains all the information needed to pick up
values in the source CPU_grid an affect them in the destination
CPU_grid. It handles all the hydrodynamics variables at once, in order
to lower the communication cost.
 */

struct jcommunicator {
  int tag;			/**< Tag used by MPI. Necessarily
				   'int' and not 'long' as this has to
				   correspond to MPI implementation */
  long imin_dest[3]; /**< lower boundary of the communicator on the
				     destination CPU_Grid. The lower
				     boundary is included.
		     */
  long imax_dest[3]; /**< upper boundary of the communicator on the
				     destination CPU_Grid. The upper
				     boundary is excluded.
		     */
  long imin_src[3]; /**< lower boundary of the communicator on the
				     source CPU_Grid. The lower
				     boundary is included.
		    */
  long imax_src[3]; /**< upper boundary of the communicator on the
				     source CPU_Grid. The lower
				     boundary is included.
		    */
  long dest_level; /**< Refinement level of the destination grid */
  long src_level; /**< Refinement level of the source grid */
  long grid_src;  /**< Logical number of the parent of the source grid */
  long grid_dest; /**< Logical number of the parent of the destination grid */
  long facedim; /**< Dimension perpendicular to the communicator */
  long faceside; /**< INF or SUP, relative to the destination */
  long CPU_src, CPU_dest; /**< processing element number of source and destination */
  long nb_src, nb_dest; /**< Logical numbers of the source and destination grids */
  long Imin[3], Imax[3]; /**< Corners of the communicator in absolute coordinates */
  long type;			/**< Should be one of these values: GHOST, MEAN, FLUX */
  real *buffer; /**< Container for the data to be transferred */
  tGrid_CPU *srcg, *destg; /**< Pointers to the descriptors of the
			      source and destination CPU_Grids */
  long size; /**< size of communicator, for the destination mesh. The
		size of a communicator is always its size for the
		destination mesh. It is expressed in number of
		zones. */
  struct jcommunicator *next; /**< Pointer to the next communicator */
  struct jcommunicator *prev; /**< Pointer to the previous communicator */
#ifdef GPU
  struct cudaMemcpy3DParms OnSrc;
  struct cudaMemcpy3DParms OnDst;
  int dz;
  int yzsize;
  struct cudaPitchedPtr bufferGPU;
#endif
};


/** Hash tables of communicators are built to find efficiently
redundant communicators. A given hash table groups all the
communicators that have the same source CPU_Grid. The tables
are handled as chained list. */

struct commhash {
  struct jcommunicator *com; /**< Communicator of the current hash table element */
  struct commhash *next;    /**< Pointer to the next element of the hash table */
  struct commhash *prev;    /**< Pointer to the previous element of the hash table */
};

/** Structure that handles scalar hydrodynamics fields (eg density) */

struct scalarfield {
  char *Name;      /**< Name of the field (eg "density", "energy") */
  real *Field;     /**< Pointer to the array that contains the field */
  real *Field_gpu;
  boolean fresh_cpu;
  boolean fresh_gpu;
  tGrid_CPU *desc; /**< Pointer to the associated descriptor of the CPU_Grid */
  int level;
  size_t pitch;
  long   ncell[3]; 
};


/** Structure that handles vectorial hydrodynamics fields (eg velocity) */

struct vectorfield {
  char *Name;       /**< Name of the field (eg "velocity") */	       
  real *Field[3];   /**< Pointer to the array that contains the field */
  real *Field_gpu[3];
  boolean fresh_cpu[3];
  boolean fresh_gpu[3];
  tGrid_CPU *desc;  /**< Pointer to the associated descriptor of the CPU_Grid */
  int level;
  size_t pitch[3];
  long   gncell[3];
};

struct fluidpatch {
  char *Name;			/**< The name of the fluid (eg gas, or dust, etc) */
  int Fluidtype;
  ScalarField *Density;         /**< Pointer to the field \f$\rho\f$*/
  real *StartField;             /**< Pointer to the beginning of the first array (they are all contiguous in memory) */
  VectorField *Velocity;        /**< Pointer to the velocity field */
  ScalarField *Energy;          /**< Pointer to the volumic energy field */
  VectorField *V_temp;          /**< Pointer to the pivotal velocity field */        
  real *Ptr[20];  /**< Short cut to fields */
#ifdef GPU
  struct cudaPitchedPtr Ptr_gpu[20];		/**< Short cut to fields over gpu */
#endif
  real *Fluxes[3][2];           /**< Pointers to the arrays that handles all fluxes at the grid boundary */
  real **FluxesGPU;
  tGrid_CPU      *desc;		/**< Pointer to the descriptor associated to this fluid patch  */
  boolean PotentialSet;		/**< Indicates whether the potential has already been evaluated  */
  FluidPatch *next;		/**< For multifluid calculations, pointer to the next fluid  */
  FluidPatch *prev;		/**< For multifluid calculations, pointer to the previous fluid  */
  long FluidRank; 		/**< Fluid number for multifluid calculation  */
  boolean PreviousEradExists;
  ScalarField *Rho0;
  ScalarField *Energy0;
  ScalarField *Vx0;
  ScalarField *Vy0;
  ScalarField *Vz0;
};

struct gridfileinfo {
  long number;
  long linenumber;
  long nc_min[3];
  long nc_max[3];
  long size[3];
  real xmin[3], xmax[3];
  long level;
  long bc[6];
  long monoCPU;			/* Restart information */
  boolean last;
};

struct destblocks {
  long block[27][6];
  long shift[27][3];
  long nb;
};