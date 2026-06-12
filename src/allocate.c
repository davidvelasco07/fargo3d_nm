#include "fargo3d.h"
#include "J_jupiter.h"
Field *CreateFieldEmpty(char *name) {
  Field *field;
  field = (Field *) malloc(sizeof(Field));
  if (field == NULL) 
    prs_error("Insufficient memory for Field cloning");
  field->field_cpu = NULL;
  field->backup = NULL;
  field->secondary_backup = NULL;
  sprintf(field->name, "%s", name);
  field->next = ListOfGrids;     //Linkedlist
  ListOfGrids = field;
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
  sprintf(field->name, "%s", name);
  return field;
}

Field *CreateFieldAlias(char *name, Field *clone, int type) {
  Field *field;
  int i,j,k;
  
  field = (Field *) malloc(sizeof(Field));
  if (field == NULL) 
    prs_error("Insufficient memory for Field cloning");
  sprintf(field->name, "%s", name);
  field->field_cpu = clone->field_cpu; //Cloning fields
  field->backup = NULL;
  field->secondary_backup = NULL;
  field->next = ListOfGrids;     //Linkedlist
  ListOfGrids = field;
#ifdef GPU
  field->field_gpu = clone->field_gpu;
  field->gpu_pp = clone->gpu_pp;
  field->cpu_pp = clone->cpu_pp;
#endif

  *(field->fresh_cpu)     =  YES;
  *(field->fresh_gpu)     =  NO;
  
  field->type = type;
  masterprint("Grids %s and %s share their storage\n", clone->name, name);
  return field;
}

void CreateField(Field **ptr, char *name, int type, boolean sx, boolean sy, boolean sz, boolean reset) {
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
  
#ifndef GPU
  if (reset)
    memset(field->field_cpu,0,(Nx+2*NGHX)*(Ny+2*NGHY)*(Nz+2*NGHZ)*sizeof(real));
#endif

  field->desc = Current_Jupiter_Patch;
  field->level = Current_Level;

#ifdef GPU
 
  if (field->field_gpu == NULL){
    cudaMalloc((void**)&arr_gpu, Maxsize_gpu*sizeof(real));
    check_errors ("CreateField");
    //printf("Field %s has been created on the GPU with size %d\n", name, Maxsize_gpu);
    field->field_gpu     =  (real *)arr_gpu;
    cudaMemset( field->field_gpu,0, Maxsize_gpu*sizeof(real));
  }
  if(reset)
    cudaMemset( field->field_gpu,0, Maxsize_gpu*sizeof(real));
  *(field->fresh_gpu) = YES;
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
  field->desc = Current_Grid;
  
#ifdef GPU
  size = (field->desc->Pitch2D)*size2;
  if (field->field_gpu == NULL || size > field->size || field->desc->Pitch2D == -1){
    if(cudaMallocPitch(&arr_gpu, &pitch, size1*sizeof(real), size2) == cudaSuccess){
      //printf("Field %s level %d has been created on the GPU\n", name, field->desc->level);
      //printf("Pitch = %d bytes (%d elements) size1 %d size2 %d \n", (int)pitch, (int)(pitch/sizeof(real)),size1,size2);
      field->field_gpu = (real*)arr_gpu;
      field->pitch = pitch; //number of elements
      field->desc->Pitch2D = pitch/sizeof(real);
      field->size = (field->desc->Pitch2D)*size2;
      //printf("Maxsize2D %d\n",field->size);
    }
    else{
      printf("There was an error allocating %s on the GPU.\n", field->name);
      check_errors ("CreateField2D");
      MPI_Finalize();
      exit(1);
    }
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
    field->fresh_cpu = (boolean *) malloc(sizeof(boolean));
    field->fresh_gpu = (boolean *) malloc(sizeof(boolean));
    sprintf(field->name, "%s", name);
  } else
    field = *ptr;
#ifndef GPU
  array = (int *) realloc(field->field_cpu, sizeof(int)*(Ny+2*NGHY)*(Nz+2*NGHZ));
#else
#ifndef PINNED
  array = (int *) realloc(field->field_cpu, sizeof(int)*(Ny+2*NGHY)*(Nz+2*NGHZ));
#else
  if (field->field_cpu == NULL)
    cudaMallocHost((void**)&array,sizeof(int)*(Ny+2*NGHY)*(Nz+2*NGHZ));
#endif
#endif

  if (array == NULL) 
    prs_error("Insufficient memory for FieldInt2D creation-step2.");
  field->field_cpu = array;
  field->backup = NULL;
  field->secondary_backup = NULL;
  
  //memset(field->field_cpu,0,(Nz+2*NGHZ)*(Ny+2*NGHY)*sizeof(int));
  // Reset is not necessary for FARGO field

  //  masterprint("Field2D %s has been created\n", name);
  //Now on the GPU
  
#ifdef GPU
  if (Current_Grid->level == 0) {
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
