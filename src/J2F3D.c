#include <fargo3d.h>
#include <J_jupiter.h>


void Adapt_for_JUPITER (char *filename) {
  GridFileInfo grids[MAXGRIDS];
  tGrid *grid;
  tGrid_CPU *item;
  int i,j;

  NDIM = 0;

#ifdef X
  InvCoordNb[NDIM]=0;
  NDIM++;
#endif
#ifdef Y
  InvCoordNb[NDIM]=1;
  NDIM++;
#endif
#ifdef Z
  InvCoordNb[NDIM]=2;
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

  if (NDIM < 3) {
    INHIBREFDIM3 = TRUE;
    InvCoordNb[2] = 3-InvCoordNb[0]-InvCoordNb[1];
  }
  if (NDIM < 2) INHIBREFDIM2 = TRUE;
  Refine[0] = (INHIBREFDIM1 == TRUE ? FALSE : TRUE);
  Refine[1] = (INHIBREFDIM2 == TRUE ? FALSE : TRUE);
  Refine[2] = (INHIBREFDIM3 == TRUE ? FALSE : TRUE);
  
  /* We now invert the coord. permutation, */
  /* i.e. we take its reciprocal map */
  for (i = 0; i < 3; i++) {
    j = 0;
    while (InvCoordNb[j] != i) j++;
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

  corner_min0[_X_]=XMIN;
  corner_max0[_X_]=XMAX;
  corner_min0[_Y_]=YMIN;
  corner_max0[_Y_]=YMAX;
  corner_min0[_Z_]=ZMIN;
  corner_max0[_Z_]=ZMAX;

  Periodic[_X_] = (PERIODICX ? TRUE : FALSE);
  Periodic[_Y_] = (PERIODICY ? TRUE : FALSE);
  Periodic[_Z_] = (PERIODICZ ? TRUE : FALSE);
 
  if (AddSubPatch)
    refine();

  if (!Restart) {
    printf("\n BEFORE  ScanGridFile \n");
    FARGO_SAFE(ScanGridFile (filename));
    //    MPI_Barrier (MPI_COMM_WORLD);	/* All processes will have to read the above written file */
    printf("\n AFTER ScanGridFile \n");
  }
  else {
    printf("\n BEFORE  ReadGrids \n");
    ReadGrids (NbRestart, grids);
    printf("\n AFTER  ReadGrids \n");
    ConstructGrids (grids);
  }
  //FARGO_SAFE(ScanGridFile (filename));
  grid = GridList;
  //At this point we have built the tgrids
  while (grid != NULL) {
    splitgrid (grid);
    grid = grid->next;}
  //At this point we have built the CPUgrids
  FARGO_SAFE (BuildCommunicators ());
  
  //At this point we have built the COMMUNICATORS
  //initFfromJ (grid);
  item = Grid_CPU_list;
  /* Prepare FARGO3D's primitive variables */
  Density = CreateFieldEmpty ("gasdens");
  Energy  = CreateFieldEmpty ("gasenergy");
#ifdef X
  Vx      = CreateFieldEmpty ("gasvx");
#endif
#ifdef Y
  Vy      = CreateFieldEmpty ("gasvy");
#endif
#ifdef Z
  Vz      = CreateFieldEmpty ("gasvz");
#endif
  Pot     = CreateFieldEmpty ("potential");
#ifdef RADIATIVE_TRANSFER
  Energyrad = CreateFieldEmpty ("energyrad");
  QPlus = CreateFieldEmpty ("qplus");
#endif
#ifdef LABELED 
  Label = CreateFieldEmpty ("label");
#endif
#ifdef STOCKHOLM
  Density0 = CreateFieldEmpty2D ("rho0", YZ);
  Energy0  = CreateFieldEmpty2D ("energy0", YZ);
  Vx0      = CreateFieldEmpty2D ("vx0", YZ);
  Vy0      = CreateFieldEmpty2D ("vy0", YZ);
  Vz0      = CreateFieldEmpty2D ("vz0", YZ);
#endif
  
}

/*The following function initializes all the Fargo variables from a Jupiter tgridCPU,
 * and it's fluid patch
 */
void AdaptFieldsFromJ (tGrid_CPU *grid){   
  if (grid != Current_Jupiter_Patch) {
    int di,i,j,k;
    size_t pitch;
    Current_Jupiter_Patch = grid;
    Density->field_cpu = grid->Fluid->Density->Field;
    Density->level = grid->level;
    Density->desc = grid;
#ifdef STOCKHOLM
    Density0->field_cpu = grid->Fluid->Rho0->Field;
    Density0->desc = grid;
#endif

    for(i=0;i<3;i++){
      for(j=0;j<2;j++){
	Fluxes[i][j] = grid->Fluid->Fluxes[i][j];
	LocalBC[i][j] = (int)(grid->iface[i][j]);
      }
    }
 
#ifdef X  
    Vx->field_cpu = grid->Fluid->Velocity->Field[_X_];
    Vx->level = grid->level;
    Vx->desc = grid;
#ifdef STOCKHOLM
    Vx0->field_cpu = grid->Fluid->Vx0->Field;
    Vx0->desc = grid;
#endif
    Nx = grid->ncell[_X_];		
    XMIN = grid->corner_min[_X_];
    XMAX = grid->corner_max[_X_];
#else
    Nx = 1;
#endif

#ifdef Y
    Pitch_cpu = grid->stride[_Y_];  
    Vy->field_cpu = grid->Fluid->Velocity->Field[_Y_];
    Vy->level = grid->level;
    Vy->desc = grid;
#ifdef STOCKHOLM
    Vy0->field_cpu = grid->Fluid->Vy0->Field;
    Vy0->desc = grid;
#endif
    Ny = grid->ncell[_Y_];		
    YMIN = grid->corner_min[_Y_];
    YMAX = grid->corner_max[_Y_];
#else
    Ny = 1;
    Pitch_cpu = 0; 
#endif

#ifdef Z
    Stride = Stride_cpu = grid->stride[_Z_];  
    Vz->field_cpu = grid->Fluid->Velocity->Field[_Z_];
    Vz->level = grid->level;
    Vz->desc = grid;
#ifdef STOCKHOLM
    Vz0->field_cpu = grid->Fluid->Vz0->Field;
    Vz0->desc = grid;
#endif
    Nz = grid->ncell[_Z_];		
    ZMIN = grid->corner_min[_Z_];		/**< Absolute position (min corner) */
    ZMAX = grid->corner_max[_Z_];		/**< Absolute position (max corner) */
#else
    Stride_cpu = 0;
    Nz = 1;
#endif

    Energy->field_cpu = grid->Fluid->Energy->Field;
    Energy->level = grid->level;
    Energy->desc = grid;
#ifdef STOCKHOLM
    Energy0->field_cpu = grid->Fluid->Energy0->Field;
    Energy0->desc = grid;
#endif

    Pot->field_cpu = grid->Fluid->Potential->Field;
    Pot->level = grid->level;
    Pot->desc = grid;

#ifdef LABELED
    Label->field_cpu = grid->Fluid->Label->Field;
    Label->level = grid->level;
    Label->desc = grid;
#endif
    
#ifdef GPU
    Density->field_gpu = grid->Fluid->Density->Field_gpu;
    Density->fresh_cpu = &(grid->Fluid->Density->fresh_cpu);
    Density->fresh_gpu = &(grid->Fluid->Density->fresh_gpu);
    Pitch_gpu = grid->Pitch_gpu;
    Stride_gpu = grid->Stride_gpu;
    Pitch2D = grid->Pitch2D;
#ifdef STOCKHOLM
    Density0->field_gpu = grid->Fluid->Rho0->Field_gpu;
    Density0->fresh_cpu = &(grid->Fluid->Rho0->fresh_cpu);
    Density0->fresh_gpu = &(grid->Fluid->Rho0->fresh_gpu);
#endif
#ifdef X
    Vx->field_gpu = grid->Fluid->Velocity->Field_gpu[_X_];
    Vx->fresh_cpu = &(grid->Fluid->Velocity->fresh_cpu[_X_]);
    Vx->fresh_gpu = &(grid->Fluid->Velocity->fresh_gpu[_X_]);
#ifdef STOCKHOLM
    Vx0->field_gpu = grid->Fluid->Vx0->Field_gpu;
    Vx0->fresh_cpu = &(grid->Fluid->Vx0->fresh_cpu);
    Vx0->fresh_gpu = &(grid->Fluid->Vx0->fresh_gpu);
#endif
#endif

#ifdef Y
    Vy->field_gpu = grid->Fluid->Velocity->Field_gpu[_Y_];
    Vy->fresh_cpu = &(grid->Fluid->Velocity->fresh_cpu[_Y_]);
    Vy->fresh_gpu = &(grid->Fluid->Velocity->fresh_gpu[_Y_]);
#ifdef STOCKHOLM
    Vy0->field_gpu = grid->Fluid->Vy0->Field_gpu;
    Vy0->fresh_cpu = &(grid->Fluid->Vy0->fresh_cpu);
    Vy0->fresh_gpu = &(grid->Fluid->Vy0->fresh_gpu);
#endif
#endif

#ifdef Z
    Vz->field_gpu = grid->Fluid->Velocity->Field_gpu[_Z_];
    Vz->fresh_cpu = &(grid->Fluid->Velocity->fresh_cpu[_Z_]);
    Vz->fresh_gpu = &(grid->Fluid->Velocity->fresh_gpu[_Z_]);
#ifdef STOCKHOLM
    Vz0->field_gpu = grid->Fluid->Vz0->Field_gpu;
    Vz0->fresh_cpu = &(grid->Fluid->Vz0->fresh_cpu);
    Vz0->fresh_gpu = &(grid->Fluid->Vz0->fresh_gpu);
#endif
#endif
    Energy->field_gpu = grid->Fluid->Energy->Field_gpu;
    Energy->fresh_cpu = &(grid->Fluid->Energy->fresh_cpu);
    Energy->fresh_gpu = &(grid->Fluid->Energy->fresh_gpu);
#ifdef STOCKHOLM
    Energy0->field_gpu = grid->Fluid->Energy0->Field_gpu;
    Energy0->fresh_cpu = &(grid->Fluid->Energy0->fresh_cpu);
    Energy0->fresh_gpu = &(grid->Fluid->Energy0->fresh_gpu);
#endif
    Pot->field_gpu = grid->Fluid->Potential->Field_gpu;
    Pot->fresh_cpu = &(grid->Fluid->Potential->fresh_cpu);
    Pot->fresh_gpu = &(grid->Fluid->Potential->fresh_gpu);
#ifdef LABELED
    Label->field_gpu = grid->Fluid->Label->Field_gpu;
    Label->fresh_cpu = &(grid->Fluid->Label->fresh_cpu);
    Label->fresh_gpu = &(grid->Fluid->Label->fresh_gpu);
#endif
#endif
  
#ifdef MHD
    real* bx = Bx->field_cpu;
    real* by = By->field_cpu;
    real* bz = Bz->field_cpu;
#endif
    //  var |= EMFX|EMFY|EMFZ;
#ifdef RADIATIVE_TRANSFER
    Energyrad->field_cpu = grid->Fluid->EnergyRad->Field;
    QPlus->field_cpu = grid->Fluid->Qplus->Field;
#ifdef GPU
    Energyrad->field_gpu = grid->Fluid->EnergyRad->Field_gpu;
    Energyrad->fresh_cpu = &(grid->Fluid->EnergyRad->fresh_cpu);
    Energyrad->fresh_gpu = &(grid->Fluid->EnergyRad->fresh_gpu);
    QPlus->field_gpu = grid->Fluid->Qplus->Field_gpu;
    QPlus->fresh_cpu = &(grid->Fluid->Qplus->fresh_cpu);
    QPlus->fresh_gpu = &(grid->Fluid->Qplus->fresh_gpu);
#endif
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
  
    Dx = (Xmin[NGHX+Nx]-Xmin[NGHX])/Nx;
   
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
      jstride[2] = grid->gncell[0]*grid->gncell[1];
    else
      jstride[2] = 0;

    FARGO_SAFE(CreateFields ()); //(Re)alloc work arrays
  }
}

void FARGO_for_all_patches (void(*f)()) {
  tGrid_CPU *item, *current;
  item = Grid_CPU_list;
  current = Current_Jupiter_Patch;
  do {
    if (item->cpu == CPU_Rank) {
      FARGO_SAFE(AdaptFieldsFromJ (item));
      f();
    }
    item = item->next;
  } while (item != NULL);
  if (current != NULL)
    FARGO_SAFE(AdaptFieldsFromJ (current));
}

void FARGO_for_all_patches_level (void(*f)(), int lev) {
  tGrid_CPU *item, *current;
  item = Grid_CPU_list;
  current = Current_Jupiter_Patch;
  do {
    if (item->cpu == CPU_Rank) {
      if (item->level == lev) {
        FARGO_SAFE(AdaptFieldsFromJ (item));
        f();
      }
    }
    item = item->next;
  } while (item != NULL);
  FARGO_SAFE(AdaptFieldsFromJ (current));
}

void FARGO_for_all_patches_level_intarg (void(*f)(), int lev, int arg) {
  tGrid_CPU *item, *current;
  item = Grid_CPU_list;
  current = Current_Jupiter_Patch;
  do {
    if (item->cpu == CPU_Rank) {
      if (item->level == lev) {
        FARGO_SAFE(AdaptFieldsFromJ (item));
        f(arg);
      }
    }
    item = item->next;
  } while (item != NULL);
  FARGO_SAFE(AdaptFieldsFromJ (current));
}

void For_all_patches_fullsync_level (int lev) {
  tGrid_CPU *item, *current;
  item = Grid_CPU_list;
  current = Current_Jupiter_Patch;  
  do {
    if (item->cpu == CPU_Rank) {
      if (item->level == lev) {
        FARGO_SAFE(AdaptFieldsFromJ (item));
	FillGhosts (StandardFields() | ENERGY);
      }
    }
    item = item->next;
  } while (item != NULL);
  ExecCommUp (lev,StandardFields() | ENERGY);
  FARGO_SAFE(AdaptFieldsFromJ (current));
}


