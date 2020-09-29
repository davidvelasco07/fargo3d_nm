#include "fargo3d.h"

void CheckBC () {
	jCommunicator *com;
	tGrid_CPU *cpug;
	boolean Problem = FALSE;
  	long i,j,stride[3]={0,0,0};
  	long size,sizeng,count,m,ct,ii,jj,kk;
  	char *flag=NULL;
  	long *proj=NULL,dim1=0,dim2=0,maxsize;
  	long empty[3][2*NGH];
  	ct = CoordType;
  	cpug=Grid_CPU_list;
  	while (cpug != NULL) {
    	size = sizeng = 1;
    	maxsize=0;
    	for (i = 0; i < NDIM; i++) {
      		if (cpug->gncell[i] > maxsize) maxsize = cpug->gncell[i];
      		size *= cpug->gncell[i];
      		sizeng *= cpug->ncell[i];
      		if (cpug->ncell[i] < Nghost[i])
				pWarning ("Grid %ld has a number of zones in %s = %d, i.e. lower than Nghost=%ld!!",\
		    	cpug->number,\
		  		CoordNames[ct*3+InvCoordNb[i]],\
		  		cpug->ncell[i], Nghost[i]);
      		stride[i] = cpug->stride[i];
    	}
    	proj = (long *)realloc(proj,sizeof(long)*maxsize);
    	flag = (char *)realloc(flag,sizeof(char)*size);
    	for (m = 0; m < size; m++)
      		flag[m] = (char)0;
    	com=ComListGhost;
    	for (i = 0; i < NDIM; i++) {
      		dim1 = (i == 0);
      		dim2 = 2 - (i == 2);
      		for (j = INF; j <= SUP; j++) {
				if (cpug->iface[i][j] > 0) {
	  				for (ii = 0; ii < Nghost[i]; ii++) {
	    				for (jj = 0; jj < cpug->gncell[dim1]; jj++) {
	      					for (kk = 0; kk < cpug->gncell[dim2]; kk++) {
								m = (j==INF? ii:ii+Nghost[i]+cpug->ncell[i])*stride[i]+jj*stride[dim1]+kk*stride[dim2];
								flag[m] = (char)1;
	      					}
	    				}
	  				}
				}
      		}
    	}
    	while (com != NULL) {
      		if (com->destg == cpug) {
				for (ii = com->imin_dest[0]; ii < com->imax_dest[0]; ii++) {
	  				for (jj = com->imin_dest[1]; jj < com->imax_dest[1]; jj++) {
	    				for (kk = com->imin_dest[2]; kk < com->imax_dest[2]; kk++) {
	      					m = ii*stride[0]+jj*stride[1]+kk*stride[2];
							flag[m] = (char)1;
	    				}
	  				}
				}
      		}
      		com=com->next;
    	}
    	count = 0;
    	for (m = 0; m < size; m++)
      		count += (long)flag[m];
    	if (count < size-sizeng) {
      		Problem = TRUE;
      		mastererr ("Problem with grid %ld (level %ld)\n", cpug->number, cpug->level);
      		if (cpug->Parent->linenumber >= 0)
				mastererr ("corresponding to line %ld\n", cpug->Parent->linenumber);
      		else
				mastererr ("corresponding to template\n");
      		for (i = 0; i < NDIM; i++) {
				for (j = 0; j < cpug->gncell[i]; j++)
	  				proj[j] = 0;
				for (j = 0; j < 2*NGH; j++) 
	  				empty[i][j] = FALSE;
				dim1 = (i == 0);
				dim2 = 2 - (i == 2);
				for (ii = 0; ii < cpug->gncell[i]; ii++) {
	  				for (jj = 0; jj < cpug->gncell[dim1]; jj++) {
	    				for (kk = 0; kk < cpug->gncell[dim2]; kk++) {
	      					m = ii*stride[i]+jj*stride[dim1]+kk*stride[dim2];
	      					proj[ii] += (long)flag[m];
	    				}
	  				}
				}
				for (ii = 0; ii < (cpug->gncell[i]-cpug->ncell[i])/2; ii++) {
	  				if (proj[ii] == 0) {
	    				mastererr ("Lower face perpendicular to %s is fully undetermined (layer %ld).",\
						CoordNames[ct*3+InvCoordNb[i]], ii);
	    				empty[i][ii]=TRUE;
	  				}
				}
				for (ii = 0; ii < (cpug->gncell[i]-cpug->ncell[i])/2; ii++) {
	  				if (proj[ii+Nghost[i]+cpug->ncell[i]] == 0) {
	    				mastererr ("Upper face perpendicular to %s is fully undetermined (layer %ld).",\
						CoordNames[ct*3+InvCoordNb[i]], ii);
	    				empty[i][ii+Nghost[i]]=TRUE;
	 				}
				}
      		}
      		for (i = 0; i < NDIM; i++) {
				for (ii = 0; ii < Nghost[i]; ii++) {
	  				if (empty[i][ii]) {
	    				for (jj = 0; jj < cpug->gncell[dim1]; jj++) {
	      					for (kk = 0; kk < cpug->gncell[dim2]; kk++) {
								m = ii*stride[i]+jj*stride[dim1]+kk*stride[dim2];
							flag[m]=1;
	      					}
	    				}
	  				}
	  				if (empty[i][ii+Nghost[i]]) {
	    				for (jj = 0; jj < cpug->gncell[dim1]; jj++) {
	      					for (kk = 0; kk < cpug->gncell[dim2]; kk++) {
								m = (ii+Nghost[i]+cpug->ncell[i])*stride[i]+jj*stride[dim1]+kk*stride[dim2];
								flag[m]=1;
	      					}
	    				}
	  				}
				}
      		}
      		count = 0;
      		for (m = 0; m < size; m++)
				count += (long)flag[m];
      		if (count < size-sizeng) 
				mastererr ("Partially covered ghost zones on that mesh");
    	}
    	cpug=cpug->next;
  	}
  	free(flag);
  	free(proj);
  	if (Problem && (CPU_Number > 1))
    	mastererr ("You should restart this run with one process only to get a lighter error message");
  	if (Problem) {
#ifdef WRITEGHOSTS
    	WriteDescriptor (1);
#else
    	WriteDescriptor (0);
#endif
    	mastererr ("I write grid descriptors with number 0 (excl. ghosts)");
    	mastererr ("and number 1 (incl. ghosts) to help understand the problems.");
  	}
  	if (Problem) {
    	mastererr ("Please correct the above errors");
    	prs_exit(1);
  	}
}
