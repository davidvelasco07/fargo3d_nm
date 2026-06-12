#include "fargo3d.h"
#include "J_jupiter.h"

void WriteScalarField (ScalarField *Field, char *filename) {
  FILE *hdl;
  long gncell[3], Size, i, j, k, mg, m, stride[3], exg[3], ngh[3];
  real *buffer;
  JSINPUT(Field);
  hdl = prs_opend (filename);
  for (i = 0; i < 3; i++) {
    ngh[i] = 0;
#ifdef WRITEGHOSTS
    ngh[i] = Nghost[i];
#endif
    exg[i] = Nghost[i]-ngh[i];
    //printf("exg[%ld]=%ld\n",i,exg[i]);
  }
  for (i = 0; i < 3; i++) {
    gncell[i] = Field->desc->gncell[i];
    stride[i] = Field->desc->stride[i];
    //printf("gncell[%ld]=%ld\n",i,gncell[i]);
    //printf("stride[%ld]=%ld\n",i,stride[i]);
  }
  Size = (gncell[0]-2*exg[0])*(gncell[1]-2*exg[1])*(gncell[2]-2*exg[2]);
  m = 0;
  buffer = prs_malloc (Size*sizeof(real));
  for (k = exg[2]; k < gncell[2]-exg[2]; k++) {
    for (j = exg[1]; j < gncell[1]-exg[1]; j++) {
      for (i = exg[0]; i < gncell[0]-exg[0]; i++) {
	mg = i*stride[0]+j*stride[1]+k*stride[2];
	buffer[m] = Field->Field[mg];
	m++;
      }
    }
  }
  fwrite (buffer, sizeof(real), Size, hdl);
  free (buffer);
  fclose (hdl);
}

void WriteVectorField (VectorField *Field, char *filename) {
  FILE *hdl;
  long gncell[3], stride[3], Size, d, i, j, k, m, mg, exg[3], ngh[3];
  real *buffer;
 
  hdl = prs_opend (filename);
  for (i = 0; i < 3; i++) {
    ngh[i] = 0;
#ifdef WRITEGHOSTS
    ngh[i] = Nghost[i];
#endif
    exg[i] = Nghost[i]-ngh[i];
  }
  for (i = 0; i < 3; i++) {
    gncell[i] = Field->desc->gncell[i];
    stride[i] = Field->desc->stride[i];
  }
  Size = (gncell[0]-2*exg[0])*(gncell[1]-2*exg[1])*(gncell[2]-2*exg[2]);
  buffer = prs_malloc (Size*sizeof(real));
  for (d = 0; d < NDIM; d++) {
    JVINPUT(Field,d);
    m = 0;
    for (k = exg[2]; k < gncell[2]-exg[2]; k++) {
      for (j = exg[1]; j < gncell[1]-exg[1]; j++) {
	for (i = exg[0]; i < gncell[0]-exg[0]; i++) {
	  mg = i*stride[0]+j*stride[1]+k*stride[2];
	  buffer[m] = Field->Field[d][mg];
	  m++;
	}
      }
    }
    fwrite (buffer, sizeof(real), Size, hdl);
  }
  free (buffer);
  fclose (hdl);
}

void WriteFluid (FluidPatch *patch, long number) {
  char filename[MAXLINELENGTH];
  long i;
  ScalarField *sf=NULL;
  VectorField *vf=NULL;
  fflush (stdout);
  int options = StandardFields();
  if (options & DENS){ 
      sf = patch->Density;
      sprintf (filename, "%s%s%ld_%ld_%ld.dat",\
	       patch->Name, sf->Name, number,\
	       sf->desc->number, sf->desc->level);
      WriteScalarField (sf, filename); 
  }
  if (options & ENERGY){ 
      sf = patch->Energy;
      sprintf (filename, "%s%s%ld_%ld_%ld.dat",\
	       patch->Name, sf->Name, number,\
	       sf->desc->number, sf->desc->level);
      WriteScalarField (sf, filename); 
  }
  if (options & VX || options & VY || options & VZ){ 
      vf = patch->Velocity;
      sprintf (filename, "%s%s%ld_%ld_%ld.dat",	\
	       patch->Name, vf->Name, number,\
	       vf->desc->number, vf->desc->level);
      WriteVectorField (vf, filename); 
  }
#ifdef LABELED
  if (options & LABEL){ 
      sf = patch->Label;
      sprintf (filename, "%s%s%ld_%ld_%ld.dat",\
	       patch->Name, sf->Name, number,\
	       sf->desc->number, sf->desc->level);
      WriteScalarField (sf, filename); 
  }
#endif
}

void Write (long number) {
  tGrid_CPU *grid;
  FluidPatch *Fluid;
  masterprint ("Output (dir: %s)  #%ld/%ld...", OUTPUTDIR, number, NTOT/NINTERM);
  pInfo ("Output (dir: %s) #%ld/%ld...", OUTPUTDIR, number, NTOT/NINTERM);
  WriteDescriptor (number);
  grid = Grid_CPU_list;
  while (grid != NULL) {
    if (grid->cpu == CPU_Rank) {
      Fluid = grid->Fluid;
      while (Fluid != NULL) {
	WriteFluid (Fluid, number);
	Fluid = Fluid->next;
      }
    }
    grid=grid->next;
  }
  pInfo ("Done\n");
}
