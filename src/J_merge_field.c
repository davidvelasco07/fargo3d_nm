#include "fargo3d.h"

void merge_field (char *radix, long outputnb, long ncpu, long *numbers, long *pos, long *size, long level, long cpugridnb, long nbcomp)
{
  long i, j, k, ll, c, d;
  long ii,jj,kk,lsize[3],chunk,chunkpos;
  FILE **in, *out;
  real datachunk[MAXLENGTHONEDIM];
  char filename[MAXLINELENGTH], filedest[MAXLINELENGTH], command[MAXLINELENGTH];
  in = (FILE **)prs_malloc((size_t)(ncpu * sizeof(FILE *)));
  for (i = 0; i < ncpu; i++) {	/* First we open simultaneously the source files */
    sprintf (filename, "%s%ld_%ld_%ld.dat", radix, outputnb, numbers[i], level);
    masterprint ("Processing %s\n", filename);
    in[i] = prs_openrd(filename);
  }
  sprintf (filedest, "%s", "temp.merge.file");
  out = prs_opend(filedest);
  for (d = 0; d < nbcomp; d++) {
    ll = 0;			/* We perform a loop on the 1D index l */
    do {				/* over the monoCPU grid */
      k = ll/(size[0]*size[1]);	/* l corresponds to coordinates i,j,k */
      j = (ll-k*size[0]*size[1])/size[0];
      i = ll-k*size[0]*size[1]-j*size[0];
      c = 0;			/* What CPU-grid does that correspond to ? */
      while (!((pos[0+6*c]<=i)&&(i<pos[3+6*c])&&(pos[1+6*c]<=j)&&\
	       (j<pos[4+6*c])&&(pos[2+6*c]<=k)&&(k<pos[5+6*c]))) c++;
      ii = i-pos[0+6*c];	      /* Current position on corresponding CPU-grid */
      jj = j-pos[1+6*c];
      kk = k-pos[2+6*c];
      lsize[0] = pos[3+6*c]-pos[0+6*c]; /* Size of that CPU-grid */
      lsize[1] = pos[4+6*c]-pos[1+6*c];
      lsize[2] = pos[5+6*c]-pos[2+6*c];
      chunk = pos[3+6*c]-i;      /* Maximum size of contiguous chunk on that CPU-grid */
      ll += chunk;
      if (chunk > MAXLENGTHONEDIM) 
	prs_error ("MAXLENGTHONEDIM preprocessor variable is too small. Edit source and recompile.");
      chunkpos = ii+jj*lsize[0]+kk*lsize[0]*lsize[1]+d*lsize[0]*lsize[1]*lsize[2];
      fseek (in[c], chunkpos*sizeof(real), SEEK_SET); /* We set the current position on the */
      fread (datachunk, sizeof(real), chunk, in[c]); /* input file, then read the chunk */
      fwrite (datachunk, sizeof(real), chunk, out); /* that we append to the merged file */
    } while (ll < size[0]*size[1]*size[2]);
  }
  for (i = 0; i < ncpu; i++) {	/* Finally we close all the open files and remove the source files */
    fclose (in[i]);
    sprintf (command, "/bin/rm -f %s%s%ld_%ld_%ld.dat", OutputDir, radix, outputnb, numbers[i], level);
    system (command);
  }
  fclose (out);
  sprintf (command, "/bin/mv -f %s%s %s%s%ld_%ld_%ld.dat", OutputDir,\
	   filedest, OutputDir, radix, outputnb, cpugridnb, level);
  system (command);
}
