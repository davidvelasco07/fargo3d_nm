#include "fargo3d.h"
Field *CreateFieldEmpty(char *name) {
  Field *field;
  field = (Field *) malloc(sizeof(Field));
  if (field == NULL) 
    prs_error("Insufficient memory for Field cloning");
  
  field->field_cpu = NULL;
  field->backup = NULL;
  field->secondary_backup = NULL;
  field->name = (char *) malloc(sizeof(char) * 80);
  sprintf(field->name, "%s", name);
  field->next = ListOfGrids;     //Linkedlist
  ListOfGrids = field;
  field->desc = Current_Jupiter_Patch;
  field->level = Current_Level;
  return field;
}

Field2D *CreateFieldEmpty2D(char *name, int dim) {
  Field2D *field;
  field = (Field2D *) malloc(sizeof(Field2D));
  if (field == NULL) 
    prs_error("Insufficient memory for Field cloning");
  field->field_cpu = NULL;
  field->backup = NULL;
  field->secondary_backup = NULL;
  field->kind=dim;
  field->name = (char *) malloc(sizeof(char) * 80);
  sprintf(field->name, "%s", name);
  field->desc = Current_Jupiter_Patch;
  return field;
}

Field *CreateFieldAlias(char *name, Field *clone, int type) {
  Field *field;
  //real *array;
  char *string;
  int i,j,k;
  
  field = (Field *) malloc(sizeof(Field));
  if (field == NULL) 
    prs_error("Insufficient memory for Field cloning");
  string = (char *) malloc(sizeof(char) * 80);
  if (string == NULL)
    prs_error("Insufficient memory for Field creation-step3");
  sprintf(string, "%s", name);
  field->field_cpu = clone->field_cpu; //Cloning fields
  field->backup = NULL;
  field->secondary_backup = NULL;
  field->owner = (Field **)(field->field_cpu+(Ny+2*NGHY)*(Nx+2*NGHX)*(Nz+2*NGHZ));
  field->name = string;
  field->next = ListOfGrids;     //Linkedlist
  ListOfGrids = field;
#ifdef GPU
  field->field_gpu = clone->field_gpu;
  field->gpu_pp = clone->gpu_pp;
  field->cpu_pp = clone->cpu_pp;
#endif

  *(field->fresh_cpu)     =  YES;
  for (i = 0; i < 4; i++) {
    field->fresh_inside_contour_cpu[i] = YES;
    field->fresh_outside_contour_cpu[i] = YES;
  }
  *(field->fresh_gpu)    =  NO;
  for (i = 0; i < 4; i++) {
    field->fresh_inside_contour_gpu[i] = NO;
    field->fresh_outside_contour_gpu[i] = NO;
  }
  
  field->type = type;
  masterprint("Grids %s and %s share their storage\n", clone->name, name);
  return field;
}

Fluid *CreateFluid(char *name, int fluidtype) {
  Fluid *f;
  char fieldname[80];

  char *fluidname = (char*) malloc(sizeof(char)*MAXNAMELENGTH);
  fluidname = (char *) malloc(sizeof(char) * 80);
  sprintf(fluidname, "%s", name);

  f = (Fluid *) malloc(sizeof(Fluid));
  f->name = fluidname;
  f->Fluidtype = fluidtype;
  
  sprintf(fieldname,"%s%s",name,"dens");
  f->Density = CreateFieldEmpty(fieldname);
  sprintf(fieldname,"%s%s",name,"energy");
  f->Energy  = CreateFieldEmpty(fieldname);  
  f->VxMed = NULL;
  CreateField2D (&(f->VxMed),"VxMed", YZ, 0);  
#ifdef X
  sprintf(fieldname,"%s%s",name,"vx");
  f->Vx      = CreateFieldEmpty(fieldname);
  f->Vx_temp = CreateFieldEmpty("Vx_temp");
#ifdef COLLISIONPREDICTOR
  f->Vx_half = CreateField("Vx_half", VXTEMP, 1,0,0);
#endif
#endif
#ifdef Y
  sprintf(fieldname,"%s%s",name,"vy");
  f->Vy      = CreateFieldEmpty(fieldname);
  f->Vy_temp = CreateFieldEmpty("Vy_temp");
#ifdef COLLISIONPREDICTOR
  f->Vy_half = CreateField("Vy_half", VYTEMP, 1,0,0);
#endif
#endif
#ifdef Z
  sprintf(fieldname,"%s%s",name,"vz");
  f->Vz      = CreateFieldEmpty(fieldname);
  f->Vz_temp = CreateFieldEmpty("Vz_temp");
#ifdef COLLISIONPREDICTOR
  f->Vz_half = CreateField("Vz_half", VZTEMP, 1,0,0);
#endif
#endif

#ifdef STOCKHOLM
  f->Density0 = NULL;
  CreateField2D (&(f->Density0),"rho0", YZ, 0);  
  f->Energy0 = NULL;
  CreateField2D (&(f->Energy0),"e0", YZ, 0); 
  f->Vx0 = NULL; 
  CreateField2D (&(f->Vx0),"vx0", YZ, 0);  
  f->Vy0 = NULL;
  CreateField2D (&(f->Vy0),"vy0", YZ, 0);  
  f->Vz0 = NULL;
  CreateField2D (&(f->Vz0),"vz0", YZ, 0);  
#endif

#ifdef DRAGFORCE
  real *coeffvalues;
  coeffvalues = (real*)malloc(sizeof(real)*(3));
  for(int n=0; n<3;n++) coeffvalues[n]=0.0;
  f->Coeffval = coeffvalues;
  //Coeffval[0] = Stokesnumber;                                                                                                                         //Coeffval[1] = PaeticleSize;                                                                                                                         //Coeffval[2] = RhoSolid;
#endif
  return f;
}

void CreateField(Field **ptr, char *name, int type, boolean sx, boolean sy, boolean sz) {
  /*sx = YES ==> Field is staggered in X. Useful for determining the
    domain of each field.*/
  Field *field;
  real *array;
  boolean *barray;
  void *arr_gpu;
  int i,j,k;
  size_t pitch;
  if (*ptr == NULL) {
    *ptr = field = (Field *) malloc(sizeof(Field));
    if (field == NULL) 
      prs_error("Insufficient memory for Field creation-step1.");
    field->name = (char *) malloc(sizeof(char) * 80);
    sprintf(field->name, "%s", name);
    field->field_cpu = NULL;
    field->field_gpu = NULL;
    field->fresh_cpu = (boolean *) malloc(sizeof(boolean));
    field->fresh_gpu = (boolean *) malloc(sizeof(boolean));
    field->next = ListOfGrids;     //Linkedlist
    ListOfGrids = field;
    
  } else
    field = *ptr;
  
  if (field->field_cpu == NULL){
#ifndef GPU
  array = (real *) realloc(field->field_cpu, sizeof(real)*Maxsize_cpu+sizeof(Field*));
#else
#ifndef PINNED
  array = (real *)realloc(field->field_cpu, sizeof(real)*Maxsize_cpu+sizeof(Field*));
#else
  if(field->field_cpu == NULL)
  cudaMallocHost((void**)&array,sizeof(real)*Maxsize_cpu+sizeof(Field*));
#endif
#endif
  
  if (array == NULL) 
    prs_error("Insufficient memory for Field creation-step2.");
  field->field_cpu = array;
  memset(field->field_cpu,0,(Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ)*sizeof(real));
  }
  field->backup = NULL;
  field->secondary_backup = NULL;
  field->line_origin = __LINE__;
  strncpy (field->file_origin, __FILE__, MAXLINELENGTH-1);


  field->desc = Current_Jupiter_Patch;
  field->level = Current_Level;

#ifdef GPU
 
  if (field->field_gpu == NULL){
    cudaMalloc((void**)&arr_gpu, Maxsize_gpu*sizeof(real));
    check_errors ("CreateField");
    //printf("Field %s has been created on the GPU with size %d\n", name, Maxsize_gpu);
    field->field_gpu     =  (real *)arr_gpu;
    cudaMemset( field->field_gpu,0, Maxsize_gpu*sizeof(real));
    *(field->fresh_gpu) = YES;
  }
#ifdef DEBUG
  printf("------>>>+++Maxsize of %s = %d\n",field->name,Maxsize_gpu);
#endif
#endif

 
  field->type = type;

  if (sx)
    field->x = Xmin;
  else 
    field->x = Xmed;
  if (sy)
    field->y = Ymin;
  else 
    field->y = Ymed;
  if (sz)
    field->z = Zmin;
  else 
    field->z = Zmed;
}


void CreateField2D(Field2D **ptr, char *name, int dim, boolean reset) {
  Field2D *field;
  real *array;
  void *arr_gpu;
  int i,j,k;
  size_t pitch;
  int size1, size2, size;

  if (dim == YZ) {
    size1 = Ny+2*NGHY;
    size2 = Nz+2*NGHZ;
  }

  if (dim == ZX) {
    size1 = Nx+2*NGHX;
    size2 = Nz+2*NGHZ;
  }

  if (dim == XY) {
    size1 = Nx+2*NGHX;
    size2 = Ny+2*NGHY;
  }
  if (*ptr == NULL) {
    *ptr = field = (Field2D *) malloc(sizeof(Field));
    if (field == NULL) 
      prs_error("Insufficient memory for Field2D creation-step1.");
    field->field_cpu = NULL;
    field->field_gpu = NULL;
    field->fresh_cpu = (boolean *) malloc(sizeof(boolean));
    field->fresh_gpu = (boolean *) malloc(sizeof(boolean));
    field->size = 0;
    field->name = (char *) malloc(sizeof(char) * 80);
    sprintf(field->name, "%s", name);
  } else
    field = *ptr;

  if (field->field_cpu == NULL){
#ifndef GPU
  array = (real *) realloc(field->field_cpu, sizeof(real)*Maxsize2D_cpu);
#else
#ifndef PINNED
  array = (real *) realloc(field->field_cpu,  sizeof(real)*Maxsize2D_cpu);
#else
  if (field->field_cpu == NULL)
    cudaMallocHost((void**)&array,sizeof(real)*Maxsize2D_cpu);
#endif
#endif

  if (array == NULL) 
    prs_error("Insufficient memory for Field2D creation-step2.");
  field->field_cpu = array;
  }
  if (reset)
    memset(field->field_cpu,0,size1*size2*sizeof(real));

  //  masterprint("Field2D %s has been created\n", name);
  //Now on the GPU
  field->desc = Current_Jupiter_Patch;
  
#ifdef GPU
  //size = size1*size2;
  size = (field->desc->Pitch2D)*size2;
  if (field->field_gpu == NULL || size > field->size || field->desc->Pitch2D == -1){
    if(cudaMallocPitch(&arr_gpu, &pitch, size1*sizeof(real), size2) == cudaSuccess){
      field->field_gpu = (real*)arr_gpu;
      field->pitch = pitch; //number of elements
      field->desc->Pitch2D = pitch/sizeof(real);
      field->size = size1*size2;
      field->size = (field->desc->Pitch2D)*size2;
    }
    else{
      printf("There was an error allocating %s on the GPU.\n", field->name);
      check_errors ("CreateField2D");
      MPI_Finalize();
      exit(1);
    }
  }
  else{
    if(reset)
      cudaMemset( field->field_gpu,0, field->size*sizeof(real));
  }
  if (dim == YZ) // Backward compatibility (old 2D arrays were only YZ).
    Pitch2D = field->desc->Pitch2D;
  //If the array is not in YZ, we store its pitch in a new field of the 2D structure.
#endif
  *(field->fresh_gpu)     =  NO;
  *(field->fresh_cpu)     =  YES;
  field->kind = dim;
}

void CreateFieldInt2D(FieldInt2D** ptr, char *name) {
  FieldInt2D *field;
  int *array;
  void *arr_gpu;
  int i,j,k;
  size_t pitch;

  if (*ptr == NULL) {
    *ptr = field = (FieldInt2D *) malloc(sizeof(Field));
    if (field == NULL) 
      prs_error("Insufficient memory for FieldInt2D creation-step1.");
    field->field_cpu = NULL;
    field->field_gpu = NULL;
    field->fresh_cpu = (boolean *) malloc(sizeof(boolean));
    field->fresh_gpu = (boolean *) malloc(sizeof(boolean));
    field->name = (char *) malloc(sizeof(char) * 80);
    sprintf(field->name, "%s", name);
  } else
    field = *ptr;
#ifndef GPU
  array = (int *) realloc(field->field_cpu, sizeof(int)*Maxsize2D_cpu);
#else
#ifndef PINNED
  array = (int *) realloc(field->field_cpu, sizeof(int)*Maxsize2D_cpu);
#else
  if (field->field_cpu == NULL)
    cudaMallocHost((void**)&array,sizeof(int)*Maxsize2D_cpu);
#endif
#endif

  if (array == NULL) 
    prs_error("Insufficient memory for FieldInt2D creation-step2.");
  field->field_cpu = array;
  field->backup = NULL;
  field->secondary_backup = NULL;

  //Now on the GPU
  
#ifdef GPU
  if (Current_Jupiter_Patch->level == 0) {
    if (field->field_gpu == NULL) {
      cudaMallocPitch (&arr_gpu, &pitch, (Ny+2*NGHY)*sizeof(int), Nz+2*NGHZ);
      check_errors ("CreateFieldInt2D");
      // printf("Integer field %s has been created on the GPU\n", name);
      field->field_gpu =  (int*)arr_gpu;
      Pitch_Int_gpu = pitch/sizeof(int);
    }
  }
#endif
  *(field->fresh_gpu)     =  NO;
  *(field->fresh_cpu)     =  YES;
}

