#include "jupiter.h"

static Param VariableSet[MAXVARIABLES];
static long  VariableIndex = 0;
static char CurrentReadingDirectory[MAXLINELENGTH];
extern HashParam CommandLineParams[];
extern long nbparamCL;

void
var(name, ptr, type, necessary, deflt)
char           *name;
char           *ptr;
long             type;
long             necessary;
char           *deflt;
{
  real    valuer=0.0;
  long    valuei=0;
  boolean valueb=FALSE;
  double  temp;
  char c;
  sscanf (deflt, "%lf", &temp);
  if ((type == REAL) || (type == INT))
    valuer = (real) (temp);
  if (type == INT)
    valuei = (long) valuer;
  if (type == BOOL) {
    c = toupper(deflt[0]);
    if ((c == 'Y') || (c == 'T') || (c == 'E') || (c == '1'))
      valueb = TRUE;
  }
  strcpy(VariableSet[VariableIndex].name, name);
  VariableSet[VariableIndex].variable = ptr;
  VariableSet[VariableIndex].type = type;
  VariableSet[VariableIndex].necessary = necessary;
  VariableSet[VariableIndex].read = NO;
  if (necessary == NO) {
    if (type == INT) {
      *((long *) ptr) = valuei;
    } else if (type == REAL) {
      *((real *) ptr) = valuer;
    } else if (type == STRING) {
      strcpy (ptr, deflt);
    } else if (type == BOOL) {
      *((boolean *) ptr) = valueb;
    }
  }
  VariableIndex++;
}

void ReadRedefined () {
  long i, j, found=NO;
  long            *ptri;
  real           *ptrr;
  boolean        *ptrb;
  for (j = 0; j < nbparamCL; j++) {
    for (i = 0; i < VariableIndex; i++) {
      if (strcmp (VariableSet[i].name, CommandLineParams[j].name) == 0) {
	found = YES;
	if (VariableSet[i].read == REDEFINED)
	  pWarning ("Parameter %s on command line is specified twice",\
		    CommandLineParams[j].name);
	VariableSet[i].read = REDEFINED;
	ptri = (long *) (VariableSet[i].variable);
	ptrr = (real *) (VariableSet[i].variable);
	ptrb = (boolean *) (VariableSet[i].variable);
	if (VariableSet[i].type == INT) {
	  *ptri = CommandLineParams[j].intvalue;
	} else if (VariableSet[i].type == REAL) {
	  *ptrr = CommandLineParams[j].floatvalue;
	    } else if (VariableSet[i].type == BOOL) {
	  *ptrb = CommandLineParams[j].boolvalue;
	} else if (VariableSet[i].type == STRING) {
	      strcpy (VariableSet[i].variable, CommandLineParams[j].stringvalue);
	}
      }
    }
    if (found == NO) {
      pWarning ("Parameter %s on command line is unknow", CommandLineParams[j].name);
    }
  }
}


void
ReadVariables(filename)
char *filename;
{
  char           *wk;
  long             i, found, type;

  strcpy (CurrentReadingDirectory, filename);
  wk = strrchr(CurrentReadingDirectory, '/');
  if (wk != NULL) *(wk+1)=0;
  InitVariables();
  ReadRedefined ();
  ReadVarFile (filename);
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if ((VariableSet[i].read == NO) && (VariableSet[i].necessary == NO)) {
      if (found == NO) {
	pInfo("*** Secondary variables omitted :\n");
	found = YES;
      }
      if ((type = VariableSet[i].type) == REAL)
	pInfo("\t%s -->\t Default Value = %.5g\n", VariableSet[i].name, *((real *) VariableSet[i].variable));
      if (type == INT)
	pInfo("\t%s -->\t Default Value = %ld\n", VariableSet[i].name, *((long *) VariableSet[i].variable));
      if (type == BOOL)
	pInfo("\t%s -->\t Default Value = %s\n",\
			   VariableSet[i].name, (*((long *) VariableSet[i].variable) ? "True":"False"));
      if (type == STRING)
	pInfo("\t%s -->\t Default Value = %s\n", VariableSet[i].name, VariableSet[i].variable);
    }
  }
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if ((VariableSet[i].read == NO) && (VariableSet[i].necessary == YES)) {
      if (found == NO) {
	prs_stderr("*** Fatal error : undefined mandatory variable(s)\n");
	prs_stderr("\t(ie. without default value)\n");
	found = YES;
      }
      prs_error("\t%s\n", VariableSet[i].name);
    }
  }
  if (found == YES) prs_error ("");
  SetGlobalVariables ();
}



void
ReadVarFile(filename)
char *filename;
{
  char            nm[80], s[200],stringval[81];
  char            tmp[MAXLINELENGTH], gridfilename[MAXLINELENGTH];
  char            shortname[MAXLINELENGTH], filepath[MAXLINELENGTH];
  char           *s1, *wk;
  char            testbool;
  double           temp=0.0;
  real            valuer;
  boolean         valueb;
  long            i, found, valuei, success;
  long            *ptri;
  real           *ptrr;
  boolean        *ptrb;
  FILE           *input, *gridfile;
  double          foo=0.0, *afoo;
  static          int recurs_level=0;
  MPI_Status      stat;
  afoo = &foo;
  if ((CPU_Rank)  && (recurs_level == 0)) 
    /* Impose sequential read of variable files (otherwise */
    /* crash on Compaq HP/SC with > 16-30 processors !?!...) */
    MPI_Recv(afoo, 1, MPI_DOUBLE, CPU_Rank-1, 98, MPI_COMM_WORLD, &stat);
  input = fopen(filename, "r");
  if (input == NULL) {
    filename = strcat (CurrentReadingDirectory, filename);
    input = fopen(filename, "r");
  }
  if (input == NULL)
    prs_error ("Unable to read '%s'",filename);
  prs_msg("Reading parameter file `%s'.\n", filename);
  pInfo("Reading parameter file `%s'.\n", filename);
  strcpy (CurrentReadingDirectory, filename);
  wk = strrchr(CurrentReadingDirectory, '/');
  if (wk != NULL) *(wk+1)=0;
  while (fgets(s, 199, input) != NULL) {
    success = sscanf(s, "%s ", nm);
    if ((nm[0] != '#') && (success == 1)) {	/* # begins a comment line */
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%lf", &temp);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%80s ", stringval);
      valuer = (real) temp;
      valuei = (long) valuer;
      testbool = toupper(stringval[0]);
      if ((testbool == 'Y') || (testbool == 'T') || (testbool == 'E') || (testbool == '1'))
	valueb = TRUE;		/* Y like Yes, T like True, E like Exists */
      else
	valueb = FALSE;
      for (i = 0; i < (long)strlen(nm); i++) {
	nm[i] = (char) toupper(nm[i]);
      }
      if (strcmp(nm, "INCLUDE") == 0) {
	recurs_level++;
	ReadVarFile (stringval);
	recurs_level--;
      } else {
	found = NO;
	for (i = 0; i < VariableIndex; i++) {
	  if (strcmp(nm, VariableSet[i].name) == 0) {
	    found = YES;
	    if (VariableSet[i].read == YES) {
	      pWarning ("Variable %s defined more than once.", nm);
	    }
	    if (VariableSet[i].read == REDEFINED)
	      pInfo ("Variable %s of parameter file is redefined on command line\n", nm);
	    else {
	      VariableSet[i].read = YES;
	      ptri = (long *) (VariableSet[i].variable);
	      ptrr = (real *) (VariableSet[i].variable);
	      ptrb = (boolean *) (VariableSet[i].variable);
	      if (VariableSet[i].type == INT) {
		*ptri = valuei;
	      } else if (VariableSet[i].type == REAL) {
		*ptrr = valuer;
	      } else if (VariableSet[i].type == BOOL) {
		*ptrb = valueb;
	      } else if (VariableSet[i].type == STRING) {
		strcpy (VariableSet[i].variable, stringval);
	      }
	    }
	  }
	}
	if (found == NO) {
	  pWarning("Variable %s defined but does not exist in code.", nm);
	}
      }
    } else {
      if ((strncmp(nm, "#GRIDINFO", 9)==0) && (success==1) && (recurs_level==0)) {
	/* We flush all subsequent lines to a grid file, but only CPU 0 does that */
	/* Nevertheless, other CPUs need to flush the parameter file so they need */
	/* to go through the 'while (fgets...)' loop */
	ExtractPath (filename, filepath, shortname);
	sprintf (gridfilename, "%s/ref_%s", filepath, shortname);
	if (!CPU_Rank) {
	  gridfile = fopen (gridfilename, "w");
	  if (gridfile == NULL)
	    prs_error ("Could not write grid file %s\n", gridfilename);
	}
	while (fgets(s, 199, input) != NULL)
	  if (!CPU_Rank) fprintf (gridfile, "%s", s);
	if (!CPU_Rank) fclose (gridfile);
	pInfo ("Writing grid file %s\n", gridfilename);
	EmbeddedGridFile = YES;
	strcpy (tmp, gridfilename);
	wk = strrchr(tmp, '/');
	if (wk == NULL) wk = tmp-1;
	strcpy (EmbeddedGridFileName, wk+1);
	pInfo ("(short name %s)\n", EmbeddedGridFileName);
      }
    }
  }
  fclose (input);
  if ((CPU_Rank < CPU_Number-1)  && (recurs_level == 0))
    MPI_Send(afoo, 1, MPI_DOUBLE, CPU_Rank+1, 98, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);	/* All processes will have to read the above written file */
}

void ListVariables (stream)
     FILE *stream;
{
  long i;
  long type;
  for (i = 0; i < VariableIndex; i++) {
      if ((type = VariableSet[i].type) == REAL)
	fprintf(stream, "%s\t%.15g\n", VariableSet[i].name, *((real *) VariableSet[i].variable));
      if (type == INT)
	fprintf(stream, "%s\t%ld\n", VariableSet[i].name, *((long *) VariableSet[i].variable));
      if (type == BOOL)
	fprintf(stream, "%s\t%ld\n", VariableSet[i].name, *((boolean *) VariableSet[i].variable));
      if (type == STRING)
	fprintf(stream, "%s\t%s\n", VariableSet[i].name, VariableSet[i].variable);
    }
  fflush (stream);
}

void ListVariablesIDL (stream)
     FILE *stream;
{
  long i;
  long type;
  fprintf(stream, "input_par = { $\n");
  for (i = 0; i < VariableIndex; i++) {
      if ((type = VariableSet[i].type) == REAL)
	fprintf(stream, "%s:%.15g", VariableSet[i].name, *((real *) VariableSet[i].variable));
      if (type == INT)
	fprintf(stream, "%s:%ld", VariableSet[i].name, *((long *) VariableSet[i].variable));
      if (type == BOOL)
	fprintf(stream, "%s:%ld", VariableSet[i].name, *((boolean *) VariableSet[i].variable));
      if (type == STRING)
	fprintf(stream, "%s:'%s'", VariableSet[i].name, VariableSet[i].variable);
      if (i != VariableIndex-1) fprintf(stream, ",$\n");
    }
  fprintf(stream, "}\n");
  fflush (stream);
}

