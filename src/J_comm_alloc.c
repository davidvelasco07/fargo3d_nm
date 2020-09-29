#include "fargo3d.h"

static long CommNb=0;
MPI_Request *Req=NULL;

void Comm_Alloc (jCommunicator *start) {
      jCommunicator *com;
      long size,i,nvar,sizep;
      long cpusrc, cpudest;
      com = start;
      while (com != NULL) {
            cpusrc  = com->CPU_src; /* Beware ! size of com is always size of dest... */
            cpudest = com->CPU_dest;
            if ((cpusrc == CPU_Rank) || (cpudest == CPU_Rank)) {
                  if (cpusrc != cpudest)
                        CommNb++;
                  size = sizep = 1;
                  for (i = 0; i < NDIM; i++) {
	                  size *= (com->imax_dest[i]-com->imin_dest[i]);
	                  sizep *= (com->imax_dest[i]-com->imin_dest[i])+1;
                  }
                  com->size = size;
                  nvar = 2+NDIM;
                  /* Energy needs to be allocated even in the isothermal
	            case for the first communication */
                  if (com->type == FLUX){
	                  nvar = 1+2*NDIM; // <== Density and +/- momenta fluxes
#ifdef ADIABATIC
	                  nvar++; // <== Energy
#endif
                  }
                  com->buffer = prs_malloc (nvar*size*sizeof(real));
      
#ifdef GPU
                  JCommunicatorGPU (com);
#endif
            }
            com = com->next;
      }
      Req = (MPI_Request *)realloc(Req, (size_t)((CommNb+1)*sizeof(MPI_Request))); 
}


   
