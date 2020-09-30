#include <fargo3d.h>

void Adapt_for_JUPITER(char *filename)
{
  GridFileInfo grids[MAXGRIDS];
  tGrid *grid;
  tGrid_CPU *item;
  int i, j;

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
    printf("\n BEFORE  ScanGridFile \n");
    FARGO_SAFE(ScanGridFile(filename));
    //    MPI_Barrier (MPI_COMM_WORLD);	/* All processes will have to read the above written file */
    printf("\n AFTER ScanGridFile \n");
  }
  else
  {
    printf("\n BEFORE  ReadGrids \n");
    ReadGrids(NbRestart, grids);
    printf("\n AFTER  ReadGrids \n");
    ConstructGrids(grids);
  }
  //FARGO_SAFE(ScanGridFile (filename));
  grid = GridList;
  //At this point we have built the tgrids
  while (grid != NULL)
  {
    splitgrid(grid);
    Ngrids++;
    grid = grid->next;
  }
  printf("Number of grids = %d\n", Ngrids);

  //At this point we have built the CPUgrids
  FARGO_SAFE(BuildCommunicators());
  //At this point we have built the COMMUNICATORS
  //initFfromJ (grid);
  item = Grid_CPU_list;
}

/*The following function initializes all the Fargo variables from a Jupiter tgridCPU,
 * and it's fluid patch
 */
void AdaptFieldsFromJ(tGrid_CPU *grid)
{
  if (grid != Current_Jupiter_Patch)
  {
    int di, i, j, k;
    size_t pitch;
    Current_Jupiter_Patch = grid;
    for (i = 0; i < NFluids_per_rank; i++)
      Fluids[i] = grid->Fluids[i]; //Array of fluids of this grid

    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 2; j++)
      {
        Fluxes[i][j] = grid->fluid->Fluxes[i][j];
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

#ifdef MHD
    real *bx = Bx->field_cpu;
    real *by = By->field_cpu;
    real *bz = Bz->field_cpu;
#endif
    //  var |= EMFX|EMFY|EMFZ;

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
    Current_Grid = grid;

    jstride[0] = 1;
    if (NDIM > 1)
      jstride[1] = grid->gncell[0];
    else
      jstride[1] = 0;
    if (NDIM > 2)
      jstride[2] = grid->gncell[0] * grid->gncell[1];
    else
      jstride[2] = 0;
    //This should be changed to be called just once, the arrays defined here are global arrays with a fixed size = max_size
    //and there is no need to asign a descriptor to these fields as they are only used in FARGO3D routines (not involved in NM comms)
    FARGO_SAFE(CreateFields()); //(Re)alloc work arrays

    Ncpu_x = grid->Parent->Ncpus[1];
    Ncpu_y = grid->Parent->Ncpus[2];
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
