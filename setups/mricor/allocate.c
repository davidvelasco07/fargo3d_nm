#include "fargo3d.h"

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
  field->owner = (Field **)(field->field_cpu+(Ny+2*NGHY)*Nx*(Nz+2*NGHZ));
  field->name = string;
  field->next = ListOfGrids;     //Linkedlist
  field->zmean = NULL;
  ListOfGrids = field;
#ifdef GPU
  field->field_gpu = clone->field_gpu;
  field->gpu_pp = clone->gpu_pp;
  field->cpu_pp = clone->cpu_pp;
#endif

  field->fresh_cpu     =  YES;
  for (i = 0; i < 4; i++) {
    field->fresh_inside_contour_cpu[i] = YES;
    field->fresh_outside_contour_cpu[i] = YES;
  }
  field->fresh_gpu     =  NO;
  for (i = 0; i < 4; i++) {
    field->fresh_inside_contour_gpu[i] = NO;
    field->fresh_outside_contour_gpu[i] = NO;
  }
  
  field->type = type;
  masterprint("Grids %s and %s share their storage\n", clone->name, name);
  return field;
}

Field *CreateField(char *name, int type, boolean sx, boolean sy, boolean sz) {
  /*sx = YES ==> Field is staggered in X. Useful for determining the
    domain of each field.*/

  Field *field;
  real *array;
  void *arr_gpu;
  char *string;
  int i,j,k;
  size_t pitch;

  field = (Field *) malloc(sizeof(Field));
  if (field == NULL) 
    prs_error("Insufficient memory for Field creation-step1.");

#ifndef GPU
  array = (real *) malloc(sizeof(real)*(Ny+2*NGHY)*Nx*(Nz+2*NGHZ)+sizeof(Field*));
#else
#ifndef PINNED
  array = (real *) malloc(sizeof(real)*(Ny+2*NGHY)*Nx*(Nz+2*NGHZ)+sizeof(Field*));
#else
  cudaMallocHost((void**)&array,sizeof(real)*(Ny+2*NGHY)*Nx*(Nz+2*NGHZ)+sizeof(Field*));
#endif
#endif

  if (array == NULL) 
    prs_error("Insufficient memory for Field creation-step2.");
  string = (char *) malloc(sizeof(char) * 80);
  if (string == NULL) 
    prs_error("Insufficient memory for Field creation-step3.");
  sprintf(string, "%s", name);
  field->field_cpu = array;
  field->backup = NULL;
  field->secondary_backup = NULL;
  field->name = string;
  field->owner = (Field **)(array+(Ny+2*NGHY)*Nx*(Nz+2*NGHZ));
  *(field->owner) = field;
  field->line_origin = __LINE__;
  field->zmean = NULL;
  strncpy (field->file_origin, __FILE__, MAXLINELENGTH-1);
  
  field->next = ListOfGrids;     //Linkedlist
  ListOfGrids = field;

  i = j = k = 0;

#ifdef Z
  for (k = 0; k<Nz+2*NGHZ; k++) {
#endif
#ifdef Y
    for (j = 0; j<Ny+2*NGHY; j++) {
#endif
#ifdef X
      for (i = 0; i<Nx; i++) {
#endif
	array[l] = 0.0;
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif  
  masterprint("Field %s has been created\n", name);
  //Now on the GPU
#ifdef GPU
  if (Nx > 1 ) {
    cudaMallocPitch (&arr_gpu, &pitch, Nx*sizeof(real), (Ny+2*NGHY)*(Nz+2*NGHZ));
    field->gpu_pp = make_cudaPitchedPtr (arr_gpu, pitch, Nx, Ny+2*NGHY);
  } else {
    cudaMallocPitch (&arr_gpu, &pitch, (Ny+2*NGHY)*sizeof(real), Nz+2*NGHZ);
  }
  check_errors ("CreateField");
  field->cpu_pp = make_cudaPitchedPtr (array, Nx*sizeof(real), Nx, Ny+2*NGHY);
  printf("Field %s has been created on the GPU\n", name);

  field->fresh_cpu     =  YES;
  for (i = 0; i < 4; i++) {
    field->fresh_inside_contour_cpu[i] = YES;
    field->fresh_outside_contour_cpu[i] = YES;
  }
  field->fresh_gpu     =  NO;
  for (i = 0; i < 4; i++) {
    field->fresh_inside_contour_gpu[i] = NO;
    field->fresh_outside_contour_gpu[i] = NO;
  }

  field->field_gpu     =  (real *)arr_gpu;

  //OVERWRITING A LOT OF TIMES THE SAME VARIABLES
  Pitch_gpu            =  pitch/sizeof(real);
  Stride_gpu           =  Pitch_gpu*(Ny+2*NGHY);
  if (Nx == 1) {
    Pitch_gpu = 1;
    Stride_gpu = pitch/sizeof(real);
  }
#ifdef DEBUG
  printf("------>>>Pitch of %s = %d\n",field->name,pitch);
#endif
  Host2Dev3D(field); // Do NOT remove this
#endif
  Pitch_cpu            =  Nx;
  Stride_cpu           =  Pitch_cpu*(Ny+2*NGHY);

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

  return field;
}

Field2D *CreateField2D(char *name) {
  Field2D *field;
  real *array;
  void *arr_gpu;
  char *string;
  int i,j,k;
  size_t pitch;

  field = (Field2D *) malloc(sizeof(Field));
  if (field == NULL) 
    prs_error("Insufficient memory for Field2D creation-step1.");

#ifndef GPU
  array = (real *) malloc(sizeof(real)*(Ny+2*NGHY)*(Nz+2*NGHZ));
#else
#ifndef PINNED
  array = (real *) malloc(sizeof(real)*(Ny+2*NGHY)*(Nz+2*NGHZ));
#else
  cudaMallocHost((void**)&array,sizeof(real)*(Ny+2*NGHY)*(Nz+2*NGHZ));
#endif
#endif


  if (array == NULL) 
    prs_error("Insufficient memory for Field2D creation-step2.");
  string = (char *) malloc(sizeof(char) * 80);
  if (string == NULL) 
    prs_error("Insufficient memory for Field2D creation-step3.");
  sprintf(string, "%s", name);
  field->field_cpu = array;
  field->name = string;
  
  i = j = k = 0;

#ifdef Z
  for (k = 0; k<Nz+2*NGHZ; k++) {
#endif
#ifdef Y
    for (j = 0; j<Ny+2*NGHY; j++) {
#endif
	array[l2D] = 0.0;
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif  
  masterprint("Field2D %s has been created\n", name);
  //Now on the GPU
#ifdef GPU
  if(cudaMallocPitch(&arr_gpu, &pitch, (Ny+2*NGHY)*sizeof(real), Nz+2*NGHZ) == cudaSuccess){
    printf("Field %s has created on the GPU\n", name);
    printf("Pitch = %d bytes (%d elements)\n", (int)pitch, (int)(pitch/sizeof(real)));
    field->field_gpu = (real*)arr_gpu;
  }
  else{
    printf("There was an error allocating %s on the GPU.\n", field->name);
    check_errors ("CreateField2D");
    MPI_Finalize();
    exit(1);
  }
  field->fresh_gpu     =  NO;
  Pitch2D = pitch/sizeof(real);
#endif
  field->fresh_cpu     =  YES;
  return field;
}

FieldInt2D *CreateFieldInt2D(char *name) {
  FieldInt2D *field;
  int *array;
  void *arr_gpu;
  char *string;
  int i,j,k;
  size_t pitch;

  field = (FieldInt2D *) malloc(sizeof(Field));
  if (field == NULL) 
    prs_error("Insufficient memory for FieldInt2D creation-step1.");

#ifndef GPU
  array = (int *) malloc(sizeof(int)*(Ny+2*NGHY)*(Nz+2*NGHZ));
#else
#ifndef PINNED
  array = (int *) malloc(sizeof(int)*(Ny+2*NGHY)*(Nz+2*NGHZ));
#else
  cudaMallocHost((void**)&array,sizeof(int)*(Ny+2*NGHY)*(Nz+2*NGHZ));
#endif
#endif

  if (array == NULL) 
    prs_error("Insufficient memory for FieldInt2D creation-step2.");
  string = (char *) malloc(sizeof(char) * 80);
  if (string == NULL) 
    prs_error("Insufficient memory for FieldInt2D creation-step3.");
  sprintf(string, "%s", name);
  field->field_cpu = array;
  field->backup = NULL;
  field->secondary_backup = NULL;
  field->name = string;
  
  i = j = k = 0;

#ifdef Z
  for (k = 0; k<Nz+2*NGHZ; k++) {
#endif
#ifdef Y
    for (j = 0; j<Ny+2*NGHY; j++) {
#endif
	array[l2D] = 0;
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif  
  masterprint("Field2D %s has been created\n", name);
  //Now on the GPU
#ifdef GPU
  cudaMallocPitch (&arr_gpu, &pitch, (Ny+2*NGHY)*sizeof(int), Nz+2*NGHZ);
  check_errors ("CreateFieldInt2D");
  printf("Integer field %s has been created on the GPU\n", name);
  field->field_gpu =  (int*)arr_gpu;
#endif
  Pitch_Int_gpu = pitch/sizeof(int);
  field->fresh_cpu     =  YES;
  return field;
}
