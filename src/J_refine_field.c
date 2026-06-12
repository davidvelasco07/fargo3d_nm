#include "fargo3d.h"
#include "J_jupiter.h"

void refine_field (char radix[MAXLINELENGTH],long nboutput,long nbgrid,GridFileInfo grids[MAXGRIDS],long gtr,long cpugpatch,long ncomp,tGrid_CPU *cpatch,tGrid_CPU *cpar){
  char srcfile[MAXLINELENGTH];
  char destfile[MAXLINELENGTH];
  FILE *in;
  FILE *out;
  long chunk_length, i, j, k, ll, ip, jp, kp, lp, c, nc, resolp[3], resol[3], offset[3];
  real chunk_dest[MAXLENGTHONEDIM], chunk_src[MAXLENGTHONEDIM];
  long vartype;
  long ratio, totalsize, totalsizep, size[3], sizep[3], sizeg[3], sizegp[3], lpng, lng;
  real interm;
  vartype = _other_;
  setout  (nboutput);
  sprintf (srcfile,  "%s%ld_%ld_%ld.dat", radix, nboutput, nbgrid, grids[gtr].level);
  sprintf (destfile, "%s%ld_%ld_%ld.dat", radix, nboutput, cpugpatch, grids[0].level);
  in  = prs_openrd (srcfile);
  out = prs_opend (destfile);
  for (i = 0; i < 3; i++) {
    size[i]   = grids[0].size[i];
    sizeg[i]  = size[i]+2*Nghost[i];
    sizep[i]  = grids[gtr].size[i];
    sizegp[i] = sizep[i]+2*Nghost[i];
    resolp[i] = (grids[gtr].nc_max[i]-grids[gtr].nc_min[i])/sizep[i];
    resol[i]  = (grids[0].nc_max[i]-grids[0].nc_min[i])/size[i];
    offset[i] = (grids[0].nc_min[i]-grids[gtr].nc_min[i])/resolp[i];
  }
  ratio = (Refine[0] ? 2 : 1);
  chunk_length = size[0]*resol[0]/resolp[0];
  totalsize = size[0]*size[1]*size[2];
  totalsizep= sizep[0]*sizep[1]*sizep[2];
  for (nc = 0; nc < ncomp; nc++) {
    for (lng = 0; lng < totalsize; lng+=size[0]) {
      k = lng/(size[0]*size[1]);
      j = (lng-k*size[0]*size[1])/size[0];
      i = lng-k*size[0]*size[1]-j*size[0];
      ll = (i+Nghost[0])+(j+Nghost[1])*sizeg[0]+(k+Nghost[2])*sizeg[0]*sizeg[1];
      kp = k/(Refine[2] ? 2:1)+offset[2];
      jp = j/(Refine[1] ? 2:1)+offset[1];
      ip = i/(Refine[0] ? 2:1)+offset[0];
      lpng = ip+jp*sizep[0]+kp*sizep[0]*sizep[1];
      lp = (ip+Nghost[0])+(jp+Nghost[1])*sizegp[0]+(kp+Nghost[2])*sizegp[0]*sizegp[1];
      lpng +=  totalsizep*nc;
      fseek (in, lpng*sizeof(real), SEEK_SET);
      fread (chunk_src, sizeof(real), (size_t)(chunk_length), in);
      for (c = 0; c < chunk_length*ratio; c++)
	  chunk_dest[c] = chunk_src[c/ratio];
      fwrite (chunk_dest, sizeof(real), (size_t)(chunk_length*ratio), out);
    }
  }
  fclose (in);
  fclose (out);
}
