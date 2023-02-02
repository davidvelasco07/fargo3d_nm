#include "fargo3d.h"

void GridFileError (long linenumber, char *filename, const char *template, ...)
{
  char string[MAXLINELENGTH];
  va_list ap;
  va_start (ap, template);
  vsprintf (string, template, ap);
  va_end (ap);
  mastererr ("%s in %s at line %ld.\n", string, filename, linenumber);
}

void ScanGridFile (char *filename) {
  char gridfilename[MAXLINELENGTH];
  FILE *gridfile;
  long  success;
  FILE *input, *output;
  GridFileInfo grids[MAXGRIDS];
  long line=0, nb = 0;
  double xmin[3], xmax[3], foo;	/* 'foo' is used to check trailing parameters */
#ifdef RELNESTING
  double xmed[3], xlen[3];
#ifdef X
  xmed[_X_] = 0.5*(XMAX+XMIN);
  xlen[_X_] = fabs(0.5*(XMAX-XMIN));
#endif
#ifdef Y
  xmed[_Y_] = 0.5*(YMAX+YMIN);
  xlen[_Y_] = fabs(0.5*(YMAX-YMIN));
#endif
#ifdef Z
  xmed[_Z_] = 0.5*(ZMAX+ZMIN);
  xlen[_Z_] = fabs(0.5*(ZMAX-ZMIN));
#endif
#endif
  char string[MAXLINELENGTH];
  static long number=0, expect, level, i, j;
  long linenumber=0;
  boolean NestedMeshesSetup = NO;

  /* Include template mesh as first patch */
  grids[line].number = number++;
  grids[line].linenumber = -1; //Template has no specific line number after GRIDINFO
  grids[line].xmax[_X_] = grids[line].xmin[_X_] = XMIN;
  grids[line].xmax[_Y_] = grids[line].xmin[_Y_] = YMIN;
  grids[line].xmax[_Z_] = grids[line].xmin[_Z_] = ZMIN;
#ifdef X
  grids[line].xmax[_X_] = XMAX;
#endif
#ifdef Y
  grids[line].xmax[_Y_] = YMAX;
#endif
#ifdef Z
  grids[line].xmax[_Z_] = ZMAX;
#endif

  grids[line].level = 0;
  grids[line].last = TRUE;
  if (line > 0) grids[line-1].last = FALSE;
  line++;

  input = fopen(filename, "r");

  expect = 1+2*NDIM;
  while (!feof(input)) {
    if (fgets (string, MAXLINELENGTH, input) != NULL) {
      linenumber++;
      if (strncmp(string, "#GRIDINFO", 9) == 0)
	NestedMeshesSetup = YES;
      if ((*string != '#') && (strlen(string) > 1) && NestedMeshesSetup) {
	switch (NDIM) {
	case 3: 
	  nb = sscanf (string, "%lf %lf %lf %lf %lf %lf %ld %lf", \
		       xmin, xmin+1, xmin+2, xmax, xmax+1, xmax+2, &level, &foo);
	  break;
	case 2: 
	  nb = sscanf (string, "%lf %lf %lf %lf %ld %lf", \
		       xmin, xmin+1, xmax, xmax+1, &level, &foo);
	  break;
	case 1: 
	  nb = sscanf (string, "%lf %lf %ld %lf", \
		       xmin, xmax, &level, &foo);
	  break;
	}
	if (nb < expect)
	  GridFileError (linenumber, filename,\
			 "Badly formed refinement instruction (missing values ?)");
	if (nb > expect)
	  GridFileError (linenumber, filename,\
			 "Discarding trailing refinement instructions (more than %ld values)",\
			 expect);
	if (level > MaxLevel) {
	  pInfo ("Discarding line %d in gridfile since its level (%d) is higher than %d\n",\
		 linenumber, level, MaxLevel);
	} else {
	  grids[line].number = number++;
	  grids[line].linenumber = linenumber;
	  for (i = 0; i < 3; i++)	{ /* 3, not NDIM */
#ifndef RELNESTING
	    grids[line].xmin[i] = (i < NDIM ? xmin[i]: corner_min0[i]);
	    grids[line].xmax[i] = (i < NDIM ? xmax[i]: corner_min0[i]); 
#else
	    grids[line].xmin[i] = (i < NDIM ? xmin[i]*xlen[i]+xmed[i]: corner_min0[i]);
	    grids[line].xmax[i] = (i < NDIM ? xmax[i]*xlen[i]+xmed[i]: corner_min0[i]); 
#endif 
	    /* min as well for trailing dimensions */
	  }
	  grids[line].level = level;
	  grids[line].last = TRUE;
	  if (line > 0) {
	    grids[line-1].last = FALSE;
	  }
	  line++;
	}
	if (!(line < MAXGRIDS)) {
	  mastererr ("Too many grids specified. You need to recompile the code\n");
	  mastererr  ("after increasing the value of MAXGRIDS in src/J_def.h");
	}
      }
    }
  }
  fclose (input);
  WriteDimNM(grids);
  GridAbs (grids);
  GridPos (grids);
  GridBuild (grids);
}

void WriteDimNM(GridFileInfo *grids){
  FILE *output;
  char gridfilename[MAXLINELENGTH];
  int i = 0;
  if (!CPU_Rank) {
    sprintf (gridfilename, "%s/nested_meshes.txt", OUTPUTDIR);
    output = fopen (gridfilename, "w");
    if (output == NULL)
      printf ("Could not write grid file %s\n", gridfilename);

    switch (NDIM) {
    case 1: 
      fprintf (output, "# %s_min %s_max level\n\n",	\
	       SCoordNames[CoordType*3+InvCoordNb[0]],			\
	       SCoordNames[CoordType*3+InvCoordNb[0]]);
      break;
    case 2: 
      fprintf (output, "# %s_min %s_min %s_max %s_max level\n\n", \
	       SCoordNames[CoordType*3+InvCoordNb[0]],			\
	       SCoordNames[CoordType*3+InvCoordNb[1]],			\
	       SCoordNames[CoordType*3+InvCoordNb[0]],			\
	       SCoordNames[CoordType*3+InvCoordNb[1]]);
      break;
    case 3: 
      fprintf (output, "# %s_min %s_min %s_min %s_max %s_max %s_max level\n\n", \
	       SCoordNames[CoordType*3+InvCoordNb[0]],			\
	       SCoordNames[CoordType*3+InvCoordNb[1]],			\
	       SCoordNames[CoordType*3+InvCoordNb[2]],			\
	       SCoordNames[CoordType*3+InvCoordNb[0]],			\
	       SCoordNames[CoordType*3+InvCoordNb[1]],			\
	       SCoordNames[CoordType*3+InvCoordNb[2]]);
      break;
    }
  }
  do {
    if (!CPU_Rank)
      if (NDIM == 3)
	fprintf (output, "%f %f %f %f %f %f %ld\n", grids[i].xmin[0], grids[i].xmin[1], \
		 grids[i].xmin[2], grids[i].xmax[0], grids[i].xmax[1], grids[i].xmax[2], grids[i].level);
      else if (NDIM == 2)
	fprintf (output, "%f %f %f %f %ld\n", grids[i].xmin[0], grids[i].xmin[1], \
		 grids[i].xmax[0], grids[i].xmax[1], grids[i].level);
      else
	fprintf (output, "%f %f %ld\n", grids[i].xmin[0], \
		 grids[i].xmax[0], grids[i].level);
  } while (!(grids[i++].last));
  if (!CPU_Rank) fclose (output);

}