#include <fargo3d.h>

void Adapt_for_JUPITER(char *filename)
{
  GridFileInfo grids[MAXGRIDS];
  tGrid *grid;
  tGrid_CPU *item;
  FluidPatch *fluid;
  Fluid *F3Dfluid;
  int i, j, k;

  NDIM = 0;

#ifdef X
  InvCoordNb[NDIM] = 0;
  NDIM++;
#endif
#ifdef Y
  InvCoordNb[NDIM] = 1;
  NDIM++;
#endif
#ifdef Z
  InvCoordNb[NDIM] = 2;
  NDIM++;
#endif

#ifdef SPHERICAL
  CoordType = 2;
#endif
#ifdef CYLINDRICAL
  CoordType = 1;
#endif
#ifdef CARTESIAN
  CoordType = 0;
#endif

  if (NDIM < 3)
  {
    INHIBREFDIM3 = TRUE;
    InvCoordNb[2] = 3 - InvCoordNb[0] - InvCoordNb[1];
  }
  if (NDIM < 2)
    INHIBREFDIM2 = TRUE;
  Refine[0] = (INHIBREFDIM1 == TRUE ? FALSE : TRUE);
  Refine[1] = (INHIBREFDIM2 == TRUE ? FALSE : TRUE);
  Refine[2] = (INHIBREFDIM3 == TRUE ? FALSE : TRUE);

  /* We now invert the coord. permutation, */
  /* i.e. we take its reciprocal map */
  for (i = 0; i < 3; i++)
  {
    j = 0;
    while (InvCoordNb[j] != i)
      j++;
    CoordNb[i] = j;
  }

  _X_ = _Azimuthal_ = _AZIM_ = CoordNb[0];
  _Y_ = _Radial_ = _RAD_ = CoordNb[1];
  _Z_ = _Vertical_ = _Colatitude_ = _COLAT_ = CoordNb[2];

  Ncell0[_X_] = NX;
  Ncell0[_Y_] = NY;
  Ncell0[_Z_] = NZ;

  Nghost[_X_] = NGHX;
  Nghost[_Y_] = NGHY;
  Nghost[_Z_] = NGHZ;

  corner_min0[_X_] = XMIN;
  corner_max0[_X_] = XMAX;
  corner_min0[_Y_] = YMIN;
  corner_max0[_Y_] = YMAX;
  corner_min0[_Z_] = ZMIN;
  corner_max0[_Z_] = ZMAX;

  Periodic[_X_] = (PERIODICX ? TRUE : FALSE);
  Periodic[_Y_] = (PERIODICY ? TRUE : FALSE);
  Periodic[_Z_] = (PERIODICZ ? TRUE : FALSE);

  if (AddSubPatch)
    refine();

  if (!Restart)
  {
    FARGO_SAFE(ScanGridFile(filename));
    printf("\nGrid File Scanned \n");
  }
  else
  {
    ReadGrids(NbRestart, grids);
    printf("\nRead Grids \n");
    ConstructGrids(grids);
  }
  
  grid = GridList;
  //At this point we have built the tgrids
  while (grid != NULL)
  {
    splitgrid(grid);
    Ngrids++;
    printf("Grid %d N(%d,%d,%d)\n",grid->number,grid->ncell[0],grid->ncell[1],grid->ncell[2]);
    grid = grid->next;
  }

  //At this point we have built the CPUgrids
  FARGO_SAFE(BuildCommunicators());
  //At this point we have built the COMMUNICATORS
  //initFfromJ (grid);
  item = Grid_CPU_list;
  while (item != NULL)
  {
    if (item->cpu == CPU_Rank)
    {
      FARGO_SAFE(AdaptFieldsFromJ(item));
      fluid = item->fluid;
      i=NFluids_per_rank-1;
      while (fluid != NULL){
        F3Dfluid = CreateFluid(fluid->Name, fluid->Fluidtype);
        F3Dfluid->Density->field_cpu = fluid->Density->Field;
        F3Dfluid->Energy->field_cpu = fluid->Energy->Field;
        #ifdef X
        F3Dfluid->Vx->field_cpu = fluid->Velocity->Field[_X_];
        F3Dfluid->Vx_temp->field_cpu = fluid->V_temp->Field[_X_];
        #endif
        #ifdef Y
        F3Dfluid->Vy->field_cpu = fluid->Velocity->Field[_Y_];
        F3Dfluid->Vy_temp->field_cpu = fluid->V_temp->Field[_Y_];
        #endif
        #ifdef Z
        F3Dfluid->Vz->field_cpu = fluid->Velocity->Field[_Z_];
        F3Dfluid->Vz_temp->field_cpu = fluid->V_temp->Field[_Z_];
        #endif
        for (j = 0; j < 3; j++)
          for (k = 0; k < 2; k++)
            F3Dfluid->Fluxes[j][k] = fluid->Fluxes[j][k];

        #ifdef GPU
        F3Dfluid->Density->field_gpu = fluid->Density->Field_gpu;
        F3Dfluid->Density->fresh_cpu = &(fluid->Density->fresh_cpu);
        F3Dfluid->Density->fresh_gpu = &(fluid->Density->fresh_gpu);
        #ifdef X
        F3Dfluid->Vx->field_gpu = fluid->Velocity->Field_gpu[_X_];
        F3Dfluid->Vx->fresh_cpu = &(fluid->Velocity->fresh_cpu[_X_]);
        F3Dfluid->Vx->fresh_gpu = &(fluid->Velocity->fresh_gpu[_X_]);
        F3Dfluid->Vx_temp->field_gpu = fluid->V_temp->Field_gpu[_X_];
        F3Dfluid->Vx_temp->fresh_cpu = &(fluid->V_temp->fresh_cpu[_X_]);
        F3Dfluid->Vx_temp->fresh_gpu = &(fluid->V_temp->fresh_gpu[_X_]);
        #endif

        #ifdef Y
        F3Dfluid->Vy->field_gpu = fluid->Velocity->Field_gpu[_Y_];
        F3Dfluid->Vy->fresh_cpu = &(fluid->Velocity->fresh_cpu[_Y_]);
        F3Dfluid->Vy->fresh_gpu = &(fluid->Velocity->fresh_gpu[_Y_]);
        F3Dfluid->Vy_temp->field_gpu = fluid->V_temp->Field_gpu[_Y_];
        F3Dfluid->Vy_temp->fresh_cpu = &(fluid->V_temp->fresh_cpu[_Y_]);
        F3Dfluid->Vy_temp->fresh_gpu = &(fluid->V_temp->fresh_gpu[_Y_]);
        #endif

        #ifdef Z
        F3Dfluid->Vz->field_gpu = fluid->Velocity->Field_gpu[_Z_];
        F3Dfluid->Vz->fresh_cpu = &(fluid->Velocity->fresh_cpu[_Z_]);
        F3Dfluid->Vz->fresh_gpu = &(fluid->Velocity->fresh_gpu[_Z_]);
        F3Dfluid->Vz_temp->field_gpu = fluid->V_temp->Field_gpu[_Z_];
        F3Dfluid->Vz_temp->fresh_cpu = &(fluid->V_temp->fresh_cpu[_Z_]);
        F3Dfluid->Vz_temp->fresh_gpu = &(fluid->V_temp->fresh_gpu[_Z_]);
        #endif
        F3Dfluid->Energy->field_gpu = fluid->Energy->Field_gpu;
        F3Dfluid->Energy->fresh_cpu = &(fluid->Energy->fresh_cpu);
        F3Dfluid->Energy->fresh_gpu = &(fluid->Energy->fresh_gpu);
        #endif
        item->Fluids[i--] = F3Dfluid;
        F3Dfluid->FluidRank = fluid->FluidRank;
        fluid = fluid->next;
      }
    }
    item = item->next;
  }
  item = Grid_CPU_list;
}

/*The following function initializes all the Fargo variables from a Jupiter tgridCPU,
 * and it's fluid patch
 */
void SelectGrid(tGrid_CPU *grid)
{
    int di, i, j, k;
    size_t pitch;
    Current_Jupiter_Patch = grid;
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 2; j++)
      {
        LocalBC[i][j] = (int)(grid->iface[i][j]);
      }
    }
    
#ifdef X
    Nx = grid->ncell[_X_];
    XMIN = grid->corner_min[_X_];
    XMAX = grid->corner_max[_X_];
#else
    Nx = 1;
#endif

#ifdef Y
    Pitch_cpu = grid->stride[_Y_];
    Ny = grid->ncell[_Y_];
    YMIN = grid->corner_min[_Y_];
    YMAX = grid->corner_max[_Y_];
#else
    Ny = 1;
    Pitch_cpu = 0;
#endif

#ifdef Z
    Stride = Stride_cpu = grid->stride[_Z_];
    Nz = grid->ncell[_Z_];
    ZMIN = grid->corner_min[_Z_]; /**< Absolute position (min corner) */
    ZMAX = grid->corner_max[_Z_]; /**< Absolute position (max corner) */
#else
    Stride_cpu = 0;
    Nz = 1;
#endif

#ifdef GPU
    Pitch_gpu = grid->Pitch_gpu;
    Stride_gpu = grid->Stride_gpu;
    Pitch2D = grid->Pitch2D;
#endif

    Xmin = grid->Edges[_X_];
    Ymin = grid->Edges[_Y_];
    Zmin = grid->Edges[_Z_];
    Xmed = grid->Xmed;
    Ymed = grid->Ymed;
    Zmed = grid->Zmed;
    InvDiffXmed = grid->InvDiffXmed;
    InvDiffYmed = grid->InvDiffYmed;
    InvDiffZmed = grid->InvDiffZmed;
    Sxj = grid->Sxj;
    Syj = grid->Syj;
    Szj = grid->Szj;
    Sxk = grid->Sxk;
    Syk = grid->Syk;
    Szk = grid->Szk;
    InvVj = grid->InvVj;

    Dx = (Xmin[NGHX + Nx] - Xmin[NGHX]) / Nx;

#ifdef GPU
    Xmin_d = grid->Xmin_d;
    Ymin_d = grid->Ymin_d;
    Zmin_d = grid->Zmin_d;
    Sxj_d = grid->Sxj_d;
    Syj_d = grid->Syj_d;
    Szj_d = grid->Szj_d;
    Sxk_d = grid->Sxk_d;
    Syk_d = grid->Syk_d;
    Szk_d = grid->Szk_d;
    InvVj_d = grid->InvVj_d;
#endif

    Current_Level = grid->level;

    jstride[0] = 1;
    if (NDIM > 1)
      jstride[1] = grid->gncell[0];
    else
      jstride[1] = 0;
    if (NDIM > 2)
      jstride[2] = grid->gncell[0] * grid->gncell[1];
    else
      jstride[2] = 0;
      
    Ncpu_x = grid->Parent->Ncpus[1];
    Ncpu_y = grid->Parent->Ncpus[2];

    Y0 = grid->pcorner_min[1];
    Z0 = grid->pcorner_min[2];

    J = grid->cpu % Ncpu_x;
    K = grid->cpu / Ncpu_x;
}

void AdaptFieldsFromJ(tGrid_CPU *grid)
{
  if(grid != Current_Grid){
    SelectGrid(grid);
    Current_Grid = grid;
    FARGO_SAFE(CreateFields());
    //(Re)alloc work arrays
  }
}

void FARGO_for_all_patches(void (*f)())
{
  tGrid_CPU *item, *current;
  item = Grid_CPU_list;
  current = Current_Jupiter_Patch;
  do
  {
    if (item->cpu == CPU_Rank)
    {
      FARGO_SAFE(AdaptFieldsFromJ(item));
      f();
    }
    item = item->next;
  } while (item != NULL);
  if (current != NULL)
    FARGO_SAFE(AdaptFieldsFromJ(current));
}

void FARGO_for_all_patches_level(void (*f)(), int lev)
{
  tGrid_CPU *item, *current;
  item = Grid_CPU_list;
  current = Current_Jupiter_Patch;
  do
  {
    if (item->cpu == CPU_Rank)
    {
      if (item->level == lev)
      {
        FARGO_SAFE(AdaptFieldsFromJ(item));
        f();
      }
    }
    item = item->next;
  } while (item != NULL);
  FARGO_SAFE(AdaptFieldsFromJ(current));
}

void FARGO_for_all_patches_level_intarg(void (*f)(), int lev, int arg)
{
  tGrid_CPU *item, *current;
  item = Grid_CPU_list;
  current = Current_Jupiter_Patch;
  do
  {
    if (item->cpu == CPU_Rank)
    {
      if (item->level == lev)
      {
        FARGO_SAFE(AdaptFieldsFromJ(item));
        f(arg);
      }
    }
    item = item->next;
  } while (item != NULL);
  FARGO_SAFE(AdaptFieldsFromJ(current));
}

void For_all_patches_fullsync_level(int lev)
{
  tGrid_CPU *item, *current;
  item = Grid_CPU_list;
  current = Current_Jupiter_Patch;
  do
  {
    if (item->cpu == CPU_Rank)
    {
      if (item->level == lev)
      {
        FARGO_SAFE(AdaptFieldsFromJ(item));
        FillGhosts(StandardFields() | ENERGY);
      }
    }
    item = item->next;
  } while (item != NULL);
  ExecCommUp(lev, StandardFields() | ENERGY);
  FARGO_SAFE(AdaptFieldsFromJ(current));
}
