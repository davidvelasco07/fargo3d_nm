#include "fargo3d.h"

/* The following  function  dumps an  ASCII  file  to  the disk,  that
   contains the set  of minimal information that  one needs to know in
   order to draw  a sketch of the whole  mesh architecture, either  in
   physical  space,  or on  the  absolute  cartesian mesh. The  output
   meaning conforms to the specifications */

//extern real LevelDate[REFINEUPPERLIMIT];

void WriteDescriptor (long i) 
{
  char filename[MAXLINELENGTH];
  FILE *out;
  tGrid *grid = GridList;
  tGrid_CPU *cpugrid;
  long ngrid=0, ncpu, dim, j, k, pcmin[3], pcmax[3];
  int WriteGhosts=0;
#ifdef WRITEGHOSTS
  WriteGhosts=1;
#endif
  jMakeDir (i);
  if (CPU_Rank) return;
  while (grid != NULL) {ngrid++; grid = grid->next;}
  sprintf (filename, "Descriptor%ld.dat", i);
  printf ("file: %s\n",filename);
  out = prs_opend (filename);
  fprintf (out, "%.18g %ld %.18g %.18g\n", PhysicalTime, i, 0.0, 0.0);
  fprintf (out, "%d %d\n", NDIM, CoordType);
  fprintf (out, "%ld %d %ld %d\n", ngrid, CPU_Number, (WriteGhosts ? 0 : NGH), 2);
  fprintf (out, "density\n");
  fprintf (out, "velocity\n");
  grid = GridList;
  while (grid != NULL) {
    ncpu = 0;
    cpugrid = Grid_CPU_list;
    while (cpugrid != NULL) {
      if (cpugrid->Parent == grid) ncpu++;
      cpugrid = cpugrid->next;
    }
    fprintf (out, "%ld %ld %ld\n", grid->number, grid->level, ncpu);
    if (WriteGhosts) 
      fprintf (out, "%ld %ld %ld\n", grid->gncell[0],\
	       grid->gncell[1],\
	       grid->gncell[2]);
    else
      fprintf (out, "%ld %ld %ld\n", grid->ncell[0],\
	       grid->ncell[1],\
	       grid->ncell[2]);
    for (dim = 0; dim < 3; dim++) {
      for (j = INF; j <= SUP; j++) {
	fprintf(out, "%ld ", grid->BoundaryConditions[dim][j]);
      }
    }
    fprintf (out, "\n");
    for (dim = 0; dim < 3; dim++) {
      if (Nghost[dim] < 2) {	/* We write fake Edges values (so that
				   IDL procedures do not need to test
				   dimensionality) */
	for (j = 0; j < NGH; j++)
	  fprintf (out, "%.18g ", grid->Edges[dim][0]-(real)(2-j));
      }
      for (j = 0; j < grid->ncell[dim]+1+2*Nghost[dim]; j++)
	fprintf (out, "%.18g ", grid->Edges[dim][j]);
      if (Nghost[dim] < 2) {	/* We write fake Edges values (so that
				   IDL procedures do not need to test
				   dimensionality) */
	for (j = 0; j < NGH; j++)
	  fprintf (out, "%.18g ", grid->Edges[dim][grid->ncell[dim]]+(real)(j+1));
      }
      fprintf (out, "\n");
    }
    fprintf (out, "%ld %ld %ld ", grid->ncorner_min[0], grid->ncorner_min[1], grid->ncorner_min[2]);
    fprintf (out, "%ld %ld %ld\n", grid->ncorner_max[0], grid->ncorner_max[1], grid->ncorner_max[2]);
    fprintf (out, "%ld %ld %ld ", grid->gncorner_min[0], grid->gncorner_min[1], grid->gncorner_min[2]);
    fprintf (out, "%ld %ld %ld\n", grid->gncorner_max[0], grid->gncorner_max[1], grid->gncorner_max[2]);
        cpugrid = Grid_CPU_list;
    while (cpugrid != NULL) {
      if (cpugrid->Parent == grid) {
	if (WriteGhosts)
	  fprintf (out, "%ld %ld %ld %d %ld\n", cpugrid->gncell[0],\
		   cpugrid->gncell[1],\
		   cpugrid->gncell[2],\
		   cpugrid->cpu, cpugrid->number);
	else
	  fprintf (out, "%ld %ld %ld %d %ld\n", cpugrid->ncell[0],\
		   cpugrid->ncell[1],\
		   cpugrid->ncell[2],\
		   cpugrid->cpu, cpugrid->number);
	for (k = 0; k < 3; k++) {
	  pcmin[k] = cpugrid->pcorner_min[k]; /* Parent corner min. There is no notion of refinement here */
	  pcmax[k] = cpugrid->pcorner_max[k]; /* since the parent of a CPU grid is a grid at */
	  if (WriteGhosts)	/* the same level. This is why the
				   first line works whatever the value
				   of WriteGhosts. */
	    pcmax[k] += 2*Nghost[k];
	} 
	fprintf (out, "%ld %ld %ld ", pcmin[0], pcmin[1], pcmin[2]);
	fprintf (out, "%ld %ld %ld\n", pcmax[0], pcmax[1], pcmax[2]);
	for (dim = 0; dim < 3; dim++) {
	  for (j = INF; j <= SUP; j++) {
	    fprintf(out, "%ld ", cpugrid->iface[dim][j]);
	  }
	}
	fprintf(out, "\n");
      }
      cpugrid = cpugrid->next;
    }
    grid=grid->next;
  }
  fclose(out);
}
