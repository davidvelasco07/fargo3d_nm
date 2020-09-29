#include "fargo3d.h"
#include <time.h>


void J_prs_error (const char *template, ...)
{
  va_list ap;
  if (!CPU_Rank) {
    va_start (ap, template);
    vfprintf (stderr, template, ap);
    fprintf (stderr, "\nAborted.\n");
    va_end (ap);
    fflush (stderr);
  }
  prs_exit (EXIT_FAILURE);
  }

void prs_end (const char *template, ...)
{
  va_list ap;
  pInfo (template);
  if (!CPU_Rank) {
    va_start (ap, template);
    vfprintf (stdout, template, ap);
    fprintf (stdout, "\n");
    va_end (ap);
    fflush (stderr);
  }
  prs_exit (EXIT_SUCCESS);
}

void J_prs_stderr (const char *template, ...)
{
  va_list ap;
  //  pError (template); Pb at start if no parameter file... 
  if (!CPU_Rank) {
    va_start (ap, template);
    vfprintf (stderr, template, ap);
    va_end (ap);
    fflush (stderr);
  }
}

void *prs_malloc (size_t number_of_bytes)
{
  void *ptr;
  unsigned long i;
  if (number_of_bytes <= 0) 
    return NULL;
  ptr = malloc (number_of_bytes);
  if (ptr == NULL) J_prs_error ("Not enough memory.");
  for (i = 0; i < number_of_bytes; i++)
    *((char *)ptr+i) = 0;
  return ptr;
}

real **multiple_alloc_1D (int nb,int size)
{
  long i;
  real **ptr;
  ptr = prs_malloc(nb*sizeof(real *));
  ptr[0] = prs_malloc(nb*size*sizeof(real));
  for (i = 1; i < nb; i++)
    ptr[i] = ptr[i-1]+size;
  return ptr;
}

FILE *prs_open (char *filename)
{
  char fullname[MAXLINELENGTH];
  char command[MAXLINELENGTH];
  char outputdir[MAXLINELENGTH];
  char *home;
  char rcfile[MAXLINELENGTH];
  FILE *handle, *rc;
  time_t tloc;
  static boolean AlreadyWrittenLastOut = NO;
  sprintf (fullname, "%s/%s", OUTPUTDIR, filename); 
/* The 'slash' is required because this function can be called while
   the code is reading the parameter file (e.g. to open the
   'warning.log' file) and at this point a trailing slash may has not
   been appended to OUTPUTDIR */
  handle = fopen (fullname, "w");
  if (handle == NULL) {
				/* we shall try again, after creating
				   the OUTPUTDIR directory */
    sprintf (command, "mkdir -p %s", OUTPUTDIR);
    system (command);
    handle = fopen (fullname, "w");
    if (handle == NULL)
      mastererr ("I am really unable to open %s.", fullname);
  }
  return handle;
}

void setout (long number)
{
  sprintf (OutputDir, "%soutput%05ld/", OUTPUTDIR, number);
}

void jMakeDir (long number)
{
  char command[MAXLINELENGTH];
  setout (number);
  if (!CPU_Rank) {
    sprintf (command, "mkdir -p %s", OutputDir);
    system (command);
  }
  MPI_Barrier (DomainComm);
}

FILE *prs_opend (char *filename)
{
  char fullname[MAXLINELENGTH];
  FILE *handle;
  sprintf (fullname, "%s/%s", OutputDir, filename);
  handle = fopen (fullname, "w");
  if (handle == NULL) {
    mastererr ("Unable to open %s for %s\n", fullname, filename);
    exit (1);
  }
  return handle;
}

FILE *prs_openrd (char *filename)
{
  char fullname[MAXLINELENGTH];
  FILE *handle;
  sprintf (fullname, "%s%s", OutputDir, filename);
  handle = fopen (fullname, "r");
  if (handle == NULL)
    J_prs_error ("Unable to open %s.", fullname);
  return handle;
}

inline void getgridsize (tGrid_CPU *grid, long *gncell, long *stride)
{
  int k;
  for (k = 0; k < 3; k++) {	/* do not change this 3 to a NDIM. */
    gncell[k] = grid->gncell[k];
    stride[k] = grid->stride[k];
  }
}

inline void getgridsizes (tGrid_CPU *grid, long *gncell, long *stride, long *gncells, long *strides)
{
  int k;
  for (k = 0; k < 3; k++) {	/* do not change this 3 to a NDIM. */
    gncell[k] = grid->gncell[k];
    stride[k] = grid->stride[k];
  }
}

