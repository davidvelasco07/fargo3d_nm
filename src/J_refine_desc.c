#include "fargo3d.h"

void refine () {
  long output, nbgrid, ng, gtr=0, i, j, k;	/* ng : number of grids */
  long resol[3], levelmax=0;
  tGrid *grid;
  tGrid_CPU *cgrid, *descsp=NULL, *descp=NULL;
				/* gtr : grid to refine */
  GridFileInfo grids[MAXGRIDS];
  char field_name[MAXLINELENGTH];
  real pos[6];
  sscanf (SubPatchInfo, "%ld %ld %lf %lf %lf %lf %lf %lf",\
	  &output, &nbgrid, &pos[0], &pos[1], &pos[2], &pos[3], &pos[4], &pos[5]);
  masterprint ("Refining grid #%ld\n", nbgrid);
  masterprint ("at output #%ld\n", output);
  masterprint ("Intended physical limits of subpatch:\n");
  masterprint ("(%g, %g, %g) ==> (%g, %g, %g)\n",\
	      pos[0],  pos[1],  pos[2],  pos[3],  pos[4],  pos[5]);
  ReadGrids (output, grids);
  ng=0;
  
  do {
    if (grids[ng].monoCPU == nbgrid)
      gtr = ng;
    if (grids[ng].level > levelmax) levelmax = grids[ng].level;
  } while (!grids[ng++].last);
  masterprint ("The grid to refine has level %ld\n", grids[gtr].level);
  masterprint ("Its size is %ld x %ld x %ld\n",\
		      grids[gtr].size[0], grids[gtr].size[1], grids[gtr].size[2]);
  masterprint ("Its physical limits are :\n");
  masterprint ("(%g %g %g) ==> (%g %g %g)\n",\
		      grids[gtr].xmin[0], grids[gtr].xmin[1], grids[gtr].xmin[2],\
		      grids[gtr].xmax[0], grids[gtr].xmax[1], grids[gtr].xmax[2]);
  for (i = 0; i < 3; i++) {
    if (pos[i] < grids[gtr].xmin[i]) {
      masterprint ("Flushing intended lower limit on dim %ld to %g\n", i, grids[gtr].xmin[i]);
      pos[i] = grids[gtr].xmin[i];
    }
    if (pos[i+3] > grids[gtr].xmax[i]) {
      masterprint ("Flushing intended upper limit on dim %ld to %g\n", i, grids[gtr].xmax[i]);
      pos[i+3] = grids[gtr].xmax[i];
    }
  }
  gtr++;
  for (i = ng; i > 0; i--) {
    grids[i].number = grids[i-1].number;
    grids[i].last   = grids[i-1].last;
    grids[i].level  = grids[i-1].level;
    grids[i].linenumber = grids[i-1].linenumber;
    for (j = 0; j < 3; j++) {
      grids[i].xmin[j] = grids[i-1].xmin[j];
      grids[i].xmax[j] = grids[i-1].xmax[j];
      grids[i].nc_min[j] = grids[i-1].nc_min[j];
      grids[i].nc_max[j] = grids[i-1].nc_max[j];
      grids[i].size[j] = grids[i-1].size[j];
    }
    for (j = 0; j < 6; j++)
      grids[i].bc[j] = grids[i-1].bc[j];
  }
  for (i = 1; i <= ng; i++) {
    grids[i].number++;
    grids[i].last = (i == ng ? TRUE : FALSE);
    for (j = 0 ; j < 3; j++) {
      if ((levelmax == grids[gtr].level) && (Refine[j])) {
	grids[i].nc_min[j] *= 2;
	grids[i].nc_max[j] *= 2;
      }
    }
  }
  grids[0].number = 0;
  grids[0].last = FALSE;
  grids[0].level = grids[gtr].level+1;
  grids[0].linenumber = 999;
  for (i = 0; i < 3; i++) {
    grids[0].xmin[i] = pos[i];
    grids[0].xmax[i] = pos[i+3];
    resol[i] = (grids[gtr].nc_max[i]-grids[gtr].nc_min[i])/grids[gtr].size[i];
    grids[0].nc_min[i] = (long)round((real)(grids[gtr].size[i])*\
				 (grids[0].xmin[i]-grids[gtr].xmin[i])/\
				 (grids[gtr].xmax[i]-grids[gtr].xmin[i]))*\
      resol[i]+grids[gtr].nc_min[i];
    grids[0].nc_max[i] = (long)round((real)(grids[gtr].size[i])*\
				 (grids[0].xmax[i]-grids[gtr].xmin[i])/\
				 (grids[gtr].xmax[i]-grids[gtr].xmin[i])+0.01)*\
      resol[i]+grids[gtr].nc_min[i];
    if (grids[0].nc_min[i] == grids[0].nc_max[i]) {
      if (grids[0].nc_min[i] == grids[gtr].nc_min[i])
	grids[0].nc_max[i] += resol[i];
      else
	grids[0].nc_min[i] -= resol[i];
    }
    if (Refine[i]) resol[i] /= 2;
    grids[0].size[i] = (grids[0].nc_max[i]-grids[0].nc_min[i])/resol[i];
    grids[0].bc[i] = ((i >= NDIM) || (grids[gtr].nc_min[i] == grids[0].nc_min[i]) ? grids[gtr].bc[i] : 0);
    grids[0].bc[i+3] = ((i >= NDIM) || (grids[gtr].nc_max[i] == grids[0].nc_max[i]) ? grids[gtr].bc[i+3] : 0);
  }
  for (i = 0; i < 3; i++) {
    grids[0].xmin[i] = (grids[gtr].xmax[i]-grids[gtr].xmin[i])*\
      (real)(grids[0].nc_min[i]-grids[gtr].nc_min[i])/\
      (real)(grids[gtr].nc_max[i]-grids[gtr].nc_min[i])+grids[gtr].xmin[i];
    grids[0].xmax[i] = (grids[gtr].xmax[i]-grids[gtr].xmin[i])*\
      (real)(grids[0].nc_max[i]-grids[gtr].nc_min[i])/\
      (real)(grids[gtr].nc_max[i]-grids[gtr].nc_min[i])+grids[gtr].xmin[i];
  }
  masterprint ("The adopted physical limits of subpatch are : (%g %g %g) ==> (%g %g %g)\n",\
		      grids[0].xmin[0], grids[0].xmin[1], grids[0].xmin[2],\
		      grids[0].xmax[0], grids[0].xmax[1], grids[0].xmax[2]);
  for (i = 0; i < 3; i++) {
    masterprint ("Along dim %d\n", i);
    masterprint ("The grid to refine has absolute (ghost excl.) limits %d and %d\n",\
			grids[gtr].nc_min[i], grids[gtr].nc_max[i]);
    masterprint ("The subpatch has absolute (gh. excl.) limits %d and %d\n",\
			grids[0].nc_min[i], grids[0].nc_max[i]);
  }
  ConstructGrids (grids);
  grid = GridList;
  while (grid != NULL) {splitgrid (grid); grid = grid->next;}
  WriteDescriptor (output);
  cgrid = Grid_CPU_list;
  while (cgrid) {
    if (cgrid->number == nbgrid) {
      descp = cgrid;
    }
    if (cgrid->number == ng) {
      descsp = cgrid;
    }
    cgrid = cgrid->next;
  }
  
  for (k = 0; k < NbFluids; k++) {
    sprintf (field_name, "%s%s","gas", "density");
    refine_field (field_name,  output, nbgrid, grids, gtr, ng, 1L, descsp, descp);
    sprintf (field_name, "%s%s","gas", "energy");
    refine_field (field_name,  output, nbgrid, grids, gtr, ng, 1L, descsp, descp);
    sprintf (field_name, "%s%s","gas", "velocity");
    refine_field (field_name,  output, nbgrid, grids, gtr, ng, 3L, descsp, descp);
  }
  /* nbgrid : CPUgrid number of grid to refine */
				/* grids : the grid array */
				/* gtr : number of grid to refine in this array */
				/* The refined grid has number 0 */
				/* ng is the CPU-grid number of grid to refine */
  prs_end ("Subpatch creation completed.");
}
