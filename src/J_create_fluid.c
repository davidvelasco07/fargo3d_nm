#include "fargo3d.h"

ScalarField *CreateScalarField2D(tGrid_CPU *desc, char *name, real *ptr)
{
  char *Name;
  char string[MAXLINELENGTH];
  ScalarField *sf;
  long size[3], i;
  real *Field;
  void *ptr_gpu;
  size_t pitch;
  sf = prs_malloc(sizeof(ScalarField));
  for (i = 0; i < 3; i++)
    size[i] = (i > 0 ? desc->gncell[i] : 1);
  if (strlen(name) > MAXNAMELENGTH)
  {
    sprintf(string, "%s is longer than the %d characters limit.", name, MAXNAMELENGTH);
    mastererr(string);
  }
#ifdef GPU
  cudaMallocPitch(&ptr_gpu, &pitch, size[1] * sizeof(real), size[2]);
  check_errors("CreateScalarField2D");
#ifdef DEBUG
  printf("Field %s has been created on the GPU\n", name);
  printf("Pitch = %d bytes (%d elements)\n", (int)pitch, (int)(pitch / sizeof(real)));
#endif
  sf->Field_gpu = (real *)ptr_gpu;
  desc->Pitch2D = pitch / sizeof(real);
#endif
  Name = prs_malloc(sizeof(char) * (MAXNAMELENGTH + 1));
  strcpy(Name, name);
  sf->Name = Name;
  sf->desc = desc;
  Field = ptr;
  sf->Field = Field;
  sf->fresh_cpu = YES;
  sf->fresh_gpu = NO;
  return sf;
}

ScalarField *CreateScalarField(tGrid_CPU *desc, char *name, real *ptr)
{
  char *Name;
  char string[MAXLINELENGTH];
  ScalarField *sf;
  long size[3], i;
  real *Field;
  size_t pitch;
  void *ptr_gpu;
  sf = prs_malloc(sizeof(ScalarField));
  if (strlen(name) > MAXNAMELENGTH)
  {
    sprintf(string, "%s is longer than the %d characters limit.", name, MAXNAMELENGTH);
    mastererr(string);
  }
  Name = prs_malloc(sizeof(char) * (MAXNAMELENGTH + 1));
  strcpy(Name, name);
  sf->Name = Name;
  sf->desc = desc;
  for (i = 0; i < 3; i++)
    size[i] = desc->gncell[i];
  Field = ptr;
  sf->Field = Field;
  sf->level = desc->level;
#ifdef GPU
  cudaMallocPitch(&ptr_gpu, &pitch, size[0] * sizeof(real), size[1] * size[2]);
  check_errors("CreateScalarField");
#ifdef DEBUG
  printf("Field %s has been created on the GPU\n", name);
  printf("-0-0-0- Pitch = %d bytes (%d elements)\n", (int)pitch, (int)(pitch / sizeof(real)));
#endif  
  sf->Field_gpu = (real *)ptr_gpu;
  sf->pitch = pitch;
#endif
  return sf;
}

void FreeScalarField(ScalarField *sf)
{
  free(sf->Name);
  free(sf);
}

VectorField *CreateVectorField(tGrid_CPU *desc, char *name, real *ptr)
{
  char *Name;
  char string[MAXLINELENGTH];
  VectorField *vf;
  long size[3], i;
  real *Field;
  size_t pitch;
  void *ptr_gpu;
  vf = prs_malloc(sizeof(VectorField));
  if (strlen(name) > MAXNAMELENGTH)
  {
    sprintf(string, "%s is longer than the %d characters limit.", name, MAXNAMELENGTH);
    mastererr(string);
  }
  Name = prs_malloc(sizeof(char) * (MAXNAMELENGTH + 1));
  strcpy(Name, name);
  vf->Name = Name;
  vf->desc = desc;
  vf->level = desc->level;
  for (i = 0; i < 3; i++)
    size[i] = desc->gncell[i];

  Field = ptr;
#ifdef GPU
  cudaMallocPitch(&ptr_gpu, &pitch, size[0] * sizeof(real), size[1] * size[2] * NDIM);
  check_errors("CreateVectorField");
#endif
  for (i = 0; i < NDIM; i++)
  {
    vf->Field[i] = Field + i * size[0] * size[1] * size[2]; /* The velocity fields are contiguous in memory */
#ifdef GPU
    vf->Field_gpu[i] = (real *)ptr_gpu + i * pitch / sizeof(real) * size[1] * size[2];
    vf->pitch[i] = pitch;
#endif
  }
  return vf;
}

void FreeVectorField(VectorField *vf)
{
  free(vf->Name);
  free(vf);
}

FluidPatch *CreateFluidPatch(tGrid_CPU *desc, char *name, int fluidtype)
{
  char *Name;
  real *StartField;
  long dim, di, i, size[3], Size, Size2D, RealSize, dimp1, dimp2, nvar, jump;
  ScalarField *Density, *Energy;
  VectorField *Velocity, *V_temp;
  FluidPatch *patch;
  void *ptr_gpu;
  size_t pitch;
  patch = prs_malloc(sizeof(FluidPatch));
  patch->desc = desc;
  if (strlen(name) > MAXNAMELENGTH)
    mastererr("%s is longer than the %d characters limit.", name, MAXNAMELENGTH);
  Name = prs_malloc(sizeof(char) * (MAXNAMELENGTH + 1));
  strcpy(Name, name);
  patch->Name = Name;
  patch->Fluidtype = fluidtype;
  for (i = 0; i < 3; i++)
    size[i] = desc->gncell[i];
  Size = size[0] * size[1] * size[2];
  Size2D = max3(size[0] * size[1], size[0] * size[2], size[1] * size[2]);
  if (Size > Maxsize_cpu)
  {
    Maxsize_cpu = Size;
    printf("Maxsize_cpu = %d\n", Maxsize_cpu);
  }
  if (Size2D > Maxsize2D_cpu)
  {
    Maxsize2D_cpu = Size2D;
    printf("Maxsize2D_cpu = %d\n", Maxsize2D_cpu);
  }
  nvar = 2 + (2*NDIM); //Density,Energy,Velocity,V_temp
  StartField = prs_malloc(sizeof(real) * (Size * nvar));
  /* Global contiguous allocation */
  patch->StartField = StartField;
  Density = CreateScalarField(desc, "density", StartField); /* Density MUST be the first field (see function ResetPatch below) */
  Energy = CreateScalarField(desc, "energy", StartField + Size);
  Velocity = CreateVectorField(desc, "velocity", StartField + 2 * Size);
  // This leaves room to include communications for Vtemp instead of V during subcycles
  V_temp = CreateVectorField(desc, "v_temp", StartField + (2+NDIM) * Size);

  patch->Ptr[_Density_] = Density->Field;
  patch->Ptr[_Energy_] = Energy->Field;
#ifdef X
  patch->Ptr[_Vx_] = Velocity->Field[_X_];
  patch->Ptr[_Vx_temp_] = V_temp->Field[_X_];
#endif
#ifdef Y
  patch->Ptr[_Vy_] = Velocity->Field[_Y_];
  patch->Ptr[_Vy_temp_] = V_temp->Field[_Y_];
#endif
#ifdef Z
  patch->Ptr[_Vz_] = Velocity->Field[_Z_];
  patch->Ptr[_Vz_temp_] = V_temp->Field[_Z_];
#endif

#ifdef GPU
  patch->Ptr_gpu[_Density_] = make_cudaPitchedPtr(Density->Field_gpu, Density->pitch, size[0], size[1]);
  patch->Ptr_gpu[_Energy_] = make_cudaPitchedPtr(Energy->Field_gpu, Density->pitch, size[0], size[1]);
#ifdef X
  patch->Ptr_gpu[_Vx_] = make_cudaPitchedPtr(Velocity->Field_gpu[_X_], Density->pitch, size[0], size[1]);
  patch->Ptr_gpu[_Vx_temp_] = make_cudaPitchedPtr(V_temp->Field_gpu[_X_], Density->pitch, size[0], size[1]);
#endif
#ifdef Y
  patch->Ptr_gpu[_Vy_] = make_cudaPitchedPtr(Velocity->Field_gpu[_Y_], Density->pitch, size[0], size[1]);
  patch->Ptr_gpu[_Vy_temp_] = make_cudaPitchedPtr(V_temp->Field_gpu[_Y_], Density->pitch, size[0], size[1]);
#endif
#ifdef Z
  patch->Ptr_gpu[_Vz_] = make_cudaPitchedPtr(Velocity->Field_gpu[_Z_], Density->pitch, size[0], size[1]);
  patch->Ptr_gpu[_Vz_temp_] = make_cudaPitchedPtr(V_temp->Field_gpu[_Z_], Density->pitch, size[0], size[1]);
#endif

#ifdef X
  desc->Pitch_gpu = Density->pitch / sizeof(real);
  desc->Stride_gpu = desc->Pitch_gpu * (desc->gncell[_Y_]);
#else
  desc->Pitch_gpu = 1;
  desc->Stride_gpu = Density->pitch / sizeof(real);
#endif
#ifdef Y
  desc->Pitch_cpu = desc->stride[_Y_];
#else
  desc->Pitch_cpu = 0;
#endif
#ifdef Z
  desc->Stride_cpu = desc->stride[_Z_];
#else
  desc->Stride_cpu = 0;
#endif
  desc->Pitch2D = -1;

  RealSize = (Density->pitch / sizeof(real)) * size[1] * size[2];
  if (RealSize > Maxsize_gpu)
  {
    Maxsize_gpu = RealSize;
    printf("Maxsize_gpu = %d\n", Maxsize_gpu);
  }
#endif
  patch->Density = Density;
  patch->Energy = Energy;
  patch->Velocity = Velocity;
  patch->V_temp = V_temp;
#ifdef GPU
#ifdef COMMGPU
  real *Fluxes[6];
#endif
#endif

  nvar = 1 + 2 * NDIM; //Density and 2 flavors of momenta
#ifdef ADIABATIC
  nvar++; //Energy flux
#endif
  for (dim = 0; dim < 3; dim++)
  { // 3, not NDIM
    dimp1 = (dim == 0);
    dimp2 = 2 - (dim == 2);

#ifndef GPU
    patch->Fluxes[dim][INF] = prs_malloc(size[dimp1] * size[dimp2] * nvar * sizeof(real));
    patch->Fluxes[dim][SUP] = prs_malloc(size[dimp1] * size[dimp2] * nvar * sizeof(real));
#else
#ifndef COMMGPU
    patch->Fluxes[dim][INF] = prs_malloc(size[dimp1] * size[dimp2] * nvar * sizeof(real));
    patch->Fluxes[dim][SUP] = prs_malloc(size[dimp1] * size[dimp2] * nvar * sizeof(real));
#else
    cudaMalloc((void **)&patch->Fluxes[dim][INF], size[dimp1] * size[dimp2] * nvar * sizeof(real));
    cudaMalloc((void **)&patch->Fluxes[dim][SUP], size[dimp1] * size[dimp2] * nvar * sizeof(real));
    Fluxes[dim * 2 + INF] = patch->Fluxes[dim][INF];
    Fluxes[dim * 2 + SUP] = patch->Fluxes[dim][SUP];
#endif
#endif
  }

#ifdef GPU
  cudaMalloc((void **)&patch->FluxesGPU, 6 * sizeof(real));
  cudaMemcpy(patch->FluxesGPU, &Fluxes, 6 * sizeof(real), cudaMemcpyHostToDevice);
#endif
  return patch;
}

void FreeFluidPatch(FluidPatch *patch)
{
  long dim;
  free(patch->Name);
  FreeScalarField(patch->Density);
  FreeScalarField(patch->Energy);
  FreeVectorField(patch->Velocity);
  FreeVectorField(patch->V_temp);
  free(patch->StartField);
  free(patch);
}

void ResetPatch(FluidPatch *patch)
{
  long i, Size = 1;
  long nb;
#if defined(ISOTHERMAL) || defined(POLYTROPIC)
  nb = 2 + NDIM;
#else
  nb = 3 + NDIM;
#endif
  for (i = 0; i < 3; i++)
    Size *= patch->desc->gncell[i];
  for (i = 0; i < Size * nb; i++)
    patch->Density->Field[i] = 0.0;
  /* The above '2+NDIM'/'3+NDIM' is because all fields are contiguous in memory
     and 'Density' is the first one */
}
