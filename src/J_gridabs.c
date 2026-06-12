#include "fargo3d.h"
#include "J_jupiter.h"

void GridAbs (grids)
     GridFileInfo *grids;
{
  long i,le;
  i=0;
  do {
    if (grids[i].level > LevMax)
      LevMax = grids[i].level;
  } while (!(grids[i++].last));
  pInfo ("Maximum grid level is %ld\n", LevMax);
  for (i = 0; i < 3; i++) {
    le = LevMax+1;
    if (!Refine[i]) le=1;
    ncorner_min0[i] = 0;
    ncorner_max0[i] = Ncell0[i]*(1<<le);
  }
}

void GridPos (grids)
     GridFileInfo *grids;
{
  long i,j,levmax,lp,level;
  long occupyg;
  i=0;
  occupyg =  (level > 0 ? (NGH+1)/2 : 1);  //<<<<< CREO QUE ENCONTRE EL ERROR ;) ESTO HACE QUE Ncell0[0]/occuypg VALGA CERO
  do {
    level = grids[i].level;
    for (j = 0; j < 3; j++) {
      lp = (level > 0 ? level : 1);
      levmax = LevMax;
      if (!Refine[j]) { 
	levmax = 0; 
	lp = 1;
      }
      if (j < NDIM) {
	if (Refine[j]){
	  grids[i].nc_min[j] = (long)(((real)(Ncell0[j]/occupyg)*(grids[i].xmin[j]-corner_min0[j])/ \
				       (corner_max0[j]-corner_min0[j]))* \
				      (1<<(lp-1)))*(1<<(levmax+2-lp))*occupyg;
	  grids[i].nc_max[j] = (long)(((real)(Ncell0[j]/occupyg)*(grids[i].xmax[j]-corner_min0[j])/ \
				       (corner_max0[j]-corner_min0[j]))* \
				      (1<<(lp-1))+0.99999999)*(1<<(levmax+2-lp))*occupyg;
	}
	else {
	  grids[i].nc_min[j] = (long)(((real)(Ncell0[j])*(grids[i].xmin[j]-corner_min0[j])/ \
				       (corner_max0[j]-corner_min0[j]))* \
				      (1<<(lp-1)))*(1<<(levmax+2-lp))*occupyg;
	  grids[i].nc_max[j] = (long)(((real)(Ncell0[j])*(grids[i].xmax[j]-corner_min0[j])/ \
				       (corner_max0[j]-corner_min0[j]))* \
				      (1<<(lp-1))+0.99999999)*(1<<(levmax+2-lp))*occupyg;
	}
	
	if (grids[i].nc_min[j] == ncorner_min0[j])
	  grids[i].bc[j+INF*3] = 1;
	else
	  grids[i].bc[j+INF*3] = 0;
	if (grids[i].nc_max[j] == ncorner_max0[j])
	  grids[i].bc[j+SUP*3] = 1;
	else
	  grids[i].bc[j+SUP*3] = 0;
      } else {
	grids[i].nc_min[j] = 0L;
	grids[i].nc_max[j] = 2L;
	grids[i].bc[j+INF*3] = 0; //Dummy boundary conditions for
	grids[i].bc[j+SUP*3] = 0; //non-represented coordinates
      }
    }
  } while (!(grids[i++].last));
}
