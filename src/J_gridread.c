#include "fargo3d.h"

extern real LevelDate[REFINEUPPERLIMIT];

static real levdate;

void ReadGrids (long number, GridFileInfo *grids)
{
  FILE *in;
  char filename[MAXLINELENGTH];
  char line[MAXLINELENGTH];
  long foo, i, j, k, ngrid, nvar, ncpu, nb, level, ngh, size[3], bc[6];
  real rfoo;
  sprintf (filename, "%soutput%05ld/Descriptor%ld.dat", OUTPUTDIR, number, number);
  in = fopen (filename, "r");
  if (in == NULL)
    J_prs_error ("I don't find %s. I can't restart.", filename);
  printf("Reading %s\n", filename);
  //pInfo ("Reading %s\n", filename);
  fscanf (in, "%lf %ld %ld %ld", &PhysicalTime, &foo, &foo, &foo); /* Discard number: we already have it */
  fscanf (in, "%ld %ld", &foo, &foo); /* Discard NDIM & CoordType */
  fscanf (in, "%ld %ld %ld %ld", &ngrid, &ncpu, &ngh, &nvar);
  if (ncpu != 1)
    J_prs_error ("Error!! It seems that the output files are not of sequential type\n");
  for (i = 0; i < nvar; i++) 
    fscanf (in, "%s", line);	/* Discard variable names */
  for (i = 0; i < ngrid; i++) {
    fscanf (in, "%ld %ld %ld", &nb, &level, &ncpu);
    grids[nb].level = level;
    grids[nb].number = nb;
    for (j = 0; j < 3; j++)
      fscanf (in, "%ld", size+j);
    for (j = 0; j < 3; j++)
      grids[nb].size[j] = size[j];
    fscanf (in, "%ld %ld %ld %ld %ld %ld", &bc[0],&bc[3],&bc[1],&bc[4],&bc[2],&bc[5]);
    for (j = 0; j < 6; j++)
      grids[nb].bc[j] = bc[j];
    for (k = 0; k < 3; k++) {
      for (j = 0; j <= 2*ngh+size[k]; j++) {
	fscanf (in, "%lg", &rfoo);
	if (j == ngh)
	  grids[nb].xmin[k] = rfoo;
	if (j == ngh+size[k])
	  grids[nb].xmax[k] = rfoo;
      }
    }
    for (j = 0; j < 3; j++) {
      fscanf (in, "%ld", &foo);
      grids[nb].nc_min[j] = foo;
    }
    for (j = 0; j < 3; j++) {
      fscanf (in, "%ld", &foo);
      grids[nb].nc_max[j] = foo;
    }
    for (j = 0; j < 6; j++)
      fscanf (in, "%ld", &foo);	/* Discard gncorner_[min,max] information */
    for (k = 0; k < ncpu; k++) {
      for (j = 0; j < 4; j++)	
	fscanf (in, "%ld", &foo);
      fscanf (in, "%ld", &foo);
      grids[nb].monoCPU = foo;
      for (j = 0; j < 12; j++)	/* Discard CPU-grid information */
	fscanf (in, "%ld", &foo);
    }
  }
  for (i = 0; i < ngrid; i++)
    grids[i].last = FALSE;
  grids[ngrid-1].last = TRUE;
}

void ConstructGrids (GridFileInfo *grids)
{
  int i;
  GridAbs (grids);
  /* We now know LevMax */
  for (i = 0; i <= LevMax; i++)
    LevelDate[i] = levdate;
  GridBuild (grids);
}
