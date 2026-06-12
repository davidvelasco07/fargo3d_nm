#include "fargo3d.h"
#include "J_jupiter.h"

void ReadField (FluidPatch *fp,long NbRestart)
{
  char filename[MAXLINELENGTH];
  char varname[MAXLINELENGTH];
  long nvar, ndim=0, i, j, k, n, size[3], Size[3], dim, m;
  long pcmin[3], pcmax[3], gncell[3], stride[3], gsize;
  int shakehand;
  MPI_Status stat;
  FILE *in;
  real *field=NULL, *buffer;
  nvar = 3;
  getgridsize (fp->desc, gncell, stride);
  for (i = 0; i < 3; i++) {
    pcmin[i] = fp->desc->pcorner_min[i];
    pcmax[i] = fp->desc->pcorner_max[i];
    size[i] = fp->desc->ncell[i];
    Size[i] = fp->desc->Parent->ncell[i];
  }    
  gsize = gncell[0]*gncell[1]*gncell[2];
  buffer = prs_malloc(sizeof(real)*Size[0]);
  for (i = 0; i < nvar; i++) {
    switch (i) {
    case 0:
      field = fp->Density->Field;
      strcpy (varname, fp->Density->Name);
      ndim = 1;
      break;
    case 1:
      field = fp->Energy->Field;
      strcpy (varname, fp->Energy->Name);
      ndim = 1;
      break;
    case 2:
      field = fp->Velocity->Field[0];
      strcpy (varname, fp->Velocity->Name);
      ndim = NDIM;
      break;
    }
    sprintf (filename, "%soutput%05ld/%s%s%ld_%ld_%ld.dat",\
	     OUTPUTDIR, NbRestart, fp->Name, varname, NbRestart,\
	     fp->desc->Parent->monoCPU, fp->desc->level);
    //prs_msg ("Reading %s\n", filename);
    pInfo ("Reading %s\n", filename);
    if (CPU_Rank > 0) MPI_Recv (&shakehand, 1, MPI_INT, CPU_Rank-1, 43, MPI_COMM_WORLD, &stat);
    in = fopen (filename, "r");
    if (in == NULL)
      J_prs_error ("I cannot read %s. Restart impossible.", filename);
    for (dim = 0; dim < ndim; dim++) {
      for (n = 0; n < Size[1]*Size[2]; n++) {
	fread (buffer, sizeof(real), Size[0], in);
	j = n % Size[1];
	k = n / Size[1];
	if ((j >= pcmin[1]) && (k >= pcmin[2]) && (j < pcmax[1]) && (k < pcmax[2])) {
	  m = (j-pcmin[1]+Nghost[1])*stride[1]+(k-pcmin[2]+Nghost[2])*stride[2];
	  memcpy (field+m+Nghost[0]+dim*gsize, buffer+pcmin[0], size[0]*sizeof(real));
	}
      }
    }
    fclose (in);
    if (CPU_Rank < CPU_Number-1) MPI_Send (&shakehand, 1, MPI_INT, CPU_Rank+1, 43, MPI_COMM_WORLD);
  }
  free (buffer);
}
  
