#include "fargo3d.h"

real vanLeer (real slope1,real slope2)
{
  if (slope1 * slope2 <= 0.0)
    return 0;
  return (2.*slope1*slope2/(slope1+slope2));
}

int BuildFieldType(int *fieldtype, int options){
  int nb=0;
  if(options & VX)
    fieldtype[nb++] = _Vx_;
  if(options & VY) 
    fieldtype[nb++] = _Vy_;
  if(options & VZ)
    fieldtype[nb++] = _Vz_;
  if(options & DENS)
    fieldtype[nb++] = _Density_;
  if(options & ENERGY)
    fieldtype[nb++] = _Energy_;
  return nb;
}

void BuildCommSource(FluidPatch *fluid, int *fieldtype, int nvar){
	int i, field, le=0;
 	while (fluid != NULL) {
		if(fluid->FluidRank == Current_Fluid){
   			for (i = 0; i < nvar; i++) {
     			field = fieldtype[i];
     			centered[le][0] = centered[le][1] = centered[le][2] = 1;
     			if ( field == _Density_ ){
       				JSINPUT(fluid->Density);
     			}
     			else if ( field == _Energy_ ){
       				JSINPUT(fluid->Energy);
    		 	}
#ifdef X
     			else if ( field == _Vx_ ){
       				JVINPUT(fluid->Velocity,0);
       				centered[le][_X_] = 0;
     			}
#endif 
#ifdef Y
     			else if ( field == _Vy_ ){ 
       				JVINPUT(fluid->Velocity,1);
       				centered[le][_Y_] = 0;
     			}
#endif
#ifdef Z
    			else if ( field == _Vz_ ){ 
       				JVINPUT(fluid->Velocity,2);
       				centered[le][_Z_] = 0;
     			}
#endif
 				source[le++] = fluid->Ptr[field];
			}
   		}
   		fluid = fluid->next;
 	}
}

void ExecCommSame (long lev, int options){
  int fieldtype[20];
  int nvar;
  nvar = BuildFieldType(fieldtype, options);
#ifndef GPU
  ExecCommSameVar (lev, nvar, fieldtype);
#else
  ExecCommSameVar_gpu (lev, nvar, fieldtype);
#endif 
}

void ExecCommSameVar (long lev, long nvar, int *fieldtype) {	
	/* Execute intra-level communications, for level lev. These are
  	GHOST type communications */
  	jCommunicator *com;
  	long i,j,k,le,m,n,imin[3],imax[3],stride[3],d;
  	int field;
  	FluidPatch *fluid;
  	com = ComListGhost;
	while (com != NULL) {
		if ((com->dest_level == lev) && (com->src_level == lev)) {
    		if (com->CPU_src == CPU_Rank) {
				fluid = com->srcg->fluid;
				BuildCommSource(fluid,fieldtype,nvar);
				
				for (le = 0; le < 3; le++) {
	  				imin[le] = com->imin_src[le];
	  				imax[le] = com->imax_src[le];
	  				stride[le] = com->srcg->stride[le];
				}
				n=0;
				for (le=0; le < nvar; le++) {
	  				for (k = imin[2]; k < imax[2]; k++) {
	    				for (j = imin[1]; j < imax[1]; j++) {
	      					for (i = imin[0]; i < imax[0]; i++) {
								m = i*stride[0]+j*stride[1]+k*stride[2];
								com->buffer[n++] = source[le][m];
							}
	    				}
	  				}
				}
				if (com->size*nvar != n)
	  				prs_error ("Internal error: communicator size mismatch error (lev2lev)");
      		}
    	}
		com = com->next;
	}
	FARGO_SAFE(ExecComm (lev,lev,GHOST,nvar,fieldtype));
}

void ExecCommUp (long lev, int options){
  int fieldtype[20];
  int nvar;
  nvar = BuildFieldType(fieldtype, options);
#ifndef GPU
  ExecCommUpVar (lev, nvar, fieldtype);
#else
  ExecCommUpVar_gpu (lev, nvar, fieldtype);
#endif
}

void ExecCommUpVar (long lev, long nvar, int *fieldtype) /* With slope limiter */
{
  /* Execute communications from level lev to finer level lev+1. These
     are GHOST type communications. In this version we use a (multi)-linear interpolation */
  jCommunicator *com;
  long g,h,i,j,k,le,ms,mc,n,imind[3],imaxd[3],stridec[3], imaxs[3];
  long imins[3], strides[3], is[3], iid[3], size, sizec[3], msm, msp, d, isl[3];
  real s, src, srcp, srcm, loc_slope[3], ic[3];
  int field;
  FILE *hdl;
  long comp[20];
  tGrid_CPU *ds, *dd;
  long id[3], icd[3], ics[3], dn[3], side;
  long ko, jo, msl;
  real slope_i, slope_j, slope_k, value_i[9], value_j[3], frac[3];
  real vp, sr;
  FluidPatch *fluid;
  com = ComListGhost;
  while (com != NULL) {
    if ((com->dest_level == lev+1) &&		\
	(com->src_level == lev) &&		\
	(com->type == GHOST)) {
      if (com->CPU_src == CPU_Rank) {
	fluid = com->srcg->fluid;
	FARGO_SAFE(BuildCommSource(fluid,fieldtype,nvar));
	size = com->size;
	for (le = 0; le < 3; le++) {
	  imind[le] = com->imin_dest[le];
	  imaxd[le] = com->imax_dest[le];
	  imins[le] = com->imin_src[le];
	  imaxs[le] = com->imax_src[le];
	  strides[le] = com->srcg->stride[le];
	  sizec[le] = imaxd[le]-imind[le]; /* jCommunicator size */
	  stridec[le] = (le > 0 ? stridec[le-1]*sizec[le-1] : 1); /* jCommunicator stride */
	}
	n=0;
	ds = com->srcg;
	dd = com->destg;
	side = com->faceside;
	for (le=0; le<nvar; le++) {
	  for (k = imind[2]; k < imaxd[2]; k++) {
	    id[2] = k;
	    for (j = imind[1]; j < imaxd[1]; j++) {
	      id[1] = j;
	      for (i = imind[0]; i < imaxd[0]; i++) {
		id[0] = i;
		for (h = 0; h < 3; h++)
		  dn[h] = dd->dn[h];
		/* If the quantity to send is a staggered vector
		   perpendicular to the communicator and on a lower
		   face of the destination grid, we seek its position
		   on the canvas one destination cell size to the right
		   along the direction of the staggered vector, so as
		   to use the interpolated source field there (there
		   is no position mismatch between source and
		   destination, and the result at a given mesh
		   location should be independent of which kind of
		   communicator (direction X Y or Z) is used to send
		   the data). See matching comment in J_comm_exec.c */
		for (h = 0; h < 3; h++){
#if COMM == INTERNAL
		  	icd[h] = (id[h]*dn[h]+(dn[h]>>1)*centered[le][h]+(com->facedim==h?(centered[le][h]==1?0:(1-side)*dn[h]):0)+dd->gncorner_min[h]);
#endif
#if COMM == EXTERNAL
		    icd[h] = (id[h]*dn[h]+(dn[h]>>1)*centered[le][h]+(com->facedim == h ? (centered[le][h] == 1 ? 0 : side*dn[h]) : 0)+dd->gncorner_min[h]);
#endif
#if COMM == ASYMETRIC
		    icd[h] = (id[h]*dn[h]+(dn[h]>>1)*centered[le][h]+dd->gncorner_min[h]);// c like 'canvas'
#endif
		}
	    for (h = 0; h < 3; h++) {
		  	is[h] = (icd[h]-ds->gncorner_min[h])/ds->dn[h];
		  	ics[h] = (is[h]*ds->dn[h]+(ds->dn[h]>>1)*centered[le][h]+ds->gncorner_min[h]);// c like 'canvas'
		  	frac[h] = (double)(icd[h]-ics[h])/(double)(ds->dn[h]);
		}
		
		ms = 0;
		for (h = 0; h < 3; h++)
		  ms += is[h]*strides[h];
	
		for (ko=-(NDIM>2); ko<=(NDIM>2); ko++) {//-1 to +1 if 3D, 0 otherwise.
		  for (jo=-(NDIM>1); jo<=(NDIM>1); jo++) {//-1 to +1 if at least 2D, 0 otherwise.
		    msl = ms+jo*strides[1]+ko*strides[2]; //index on source mesh, local, with offset
		    src = source[le][msl];
		    srcp= source[le][msl+1];
		    srcm= source[le][msl-1];
#ifdef  COMMADAPT
		    if(centered[le][_X_] == 0){
		      vp = ds->Ymed[is[1]+jo]*sin(ds->Zmed[is[2]+ko])*OMEGAFRAME;
		      sr = sqrt(ds->Ymed[is[1]+jo]);
		      src = (src +vp)*sr;
		      srcp= (srcp+vp)*sr;
		      srcm= (srcm+vp)*sr;
		    }
#endif
		    
#ifdef NONMONOTONIC
		    if (frac[0] >= 0) {
		      value_i[jo+1+3*(ko+1)] = (1-frac[0])*src+frac[0]*srcp;
		    }
		    if (frac[0] < 0) {
		      value_i[jo+1+3*(ko+1)] = (1+frac[0])*src-frac[0]*srcm;
		    }
#else
		    slope_i = vanLeer(srcp-src,src-srcm); //Do not divide by mesh size
		    value_i[jo+1+3*(ko+1)] = src+slope_i*frac[0];
#endif
		  }
		}
		
		s = value_i[4];
		if (NDIM > 1) {
		  for (ko=-(NDIM>2); ko<=(NDIM>2); ko++) {//-1 to +1 if 3D, 0 otherwise.
		    src = value_i[1+3*(ko+1)];
		    srcp= value_i[2+3*(ko+1)];
		    srcm= value_i[3*(ko+1)];
		   
#ifdef NONMONOTONIC
		    if (frac[1] >= 0) {
		      value_j[ko+1] = (1-frac[1])*src+frac[1]*srcp;
		    }
		    if (frac[1] < 0) {
		      value_j[ko+1] = (1+frac[1])*src-frac[1]*srcm;
		    }
#else
		    slope_j = vanLeer(srcp-src,src-srcm);
		    value_j[ko+1] = src+slope_j*frac[1];
#endif
		  }
		  s = value_j[1];
		}
		if (NDIM > 2) {
		   src = value_j[1];
		   srcp= value_j[2];
		   srcm= value_j[0];
#ifdef NONMONOTONIC
		    if (frac[2] >= 0) {
		      s = (1-frac[2])*src+frac[2]*srcp;
		    }
		    if (frac[2] < 0) {
		      s = (1+frac[2])*src-frac[2]*srcm;
		    }
#else
		    slope_k = vanLeer(srcp-src,src-srcm);
		    s = src+slope_k*frac[2];
#endif
		 
		}
		mc = 0;
		for (h = 0; h < 3; h++)
		  mc += (id[h]-imind[h])*stridec[h];

		com->buffer[mc+size*le] = s;
		n++;
	      }
	    }
	  }
	}
	if (size*nvar != n)
	  prs_error ("Internal error: communicator size mismatch error");
      }
    }
    com = com->next;
  }
  FARGO_SAFE(ExecComm (lev,lev+1,GHOST,nvar,fieldtype));
}

void ExecCommDownMean (long lev, int options){
  int fieldtype[20];
  int nvar;
  nvar = BuildFieldType(fieldtype, options);
#ifndef GPU
  ExecCommDownMeanVar(lev, nvar, fieldtype);
#else
  ExecCommDownMeanVar_gpu(lev, nvar, fieldtype);
#endif
}

void ExecCommDownMeanVar (long lev, long nvar, int *fieldtype){
  /* Execute communications from level lev to coarser level lev-1, of
     type MEAN */

  jCommunicator *com;
  long i,j,k,le,ms,n,imind[3],imaxd[3], folded, ii[3], sqz[3];
  long imins[3], imaxs[3], stridec[3], strides[3],side;
  long Imin[3], Imax[3], h, dns[3], id[3], is[3], iid[3];
  real src;
  tGrid_CPU *ds, *dd;
  real coef2, coef1, coef0;
  long d,nb=2,comp[20];
  long stridecr[3], iread, iwrite;
  int field;
 
  FluidPatch *fluid;
  com = ComListMean;
 
  while (com != NULL) {
    if ((com->dest_level == lev-1) && (com->src_level == lev) && (com->type == MEAN)) {
      if (com->CPU_src == CPU_Rank) {
		fluid = com->srcg->fluid;
		BuildCommSource(fluid,fieldtype,nvar);
		ds = com->srcg;
		dd = com->destg;
		for (h = 0; h < 3; h++) {
	  		imind[h] = com->imin_dest[h];
	  		imaxd[h] = com->imax_dest[h];
	  		imins[h] = com->imin_src[h];
	  		imaxs[h] = com->imax_src[h];
	  		stridec[h] = (h > 0 ? stridec[h-1]*(imaxd[h-1]-imind[h-1]+1) : 1);
	  		stridecr[h] = (h > 0 ? stridecr[h-1]*(imaxd[h-1]-imind[h-1]) : 1);
	  		strides[h] = ds->stride[h];
	  		Imin[h] = com->Imin[h];
	  		Imax[h] = com->Imax[h];
	  		dns[h]  = ds->dn[h];
		}
		/* We first reset the buffer */
		side = com->faceside;
		for (n=0; n<com->size*nvar; n++)
	  		com->buffer[n] = 0.0;
		folded = (Refine[0] ? 2:1) * (Refine[1] ? 2:1) * (Refine[2] ? 2:1);
		for (le=0; le<nvar; le++) {
	  		for (k=0; k < (imaxd[2] - imind[2]); k += 1) {
	    		for (j=0; j < (imaxd[1] - imind[1]); j += 1) {
	      			for (i=0; i < (imaxd[0] - imind[0]); i += 1) {
						id[0] = i; id[1] = j; id[2]=k;
						n = 0;
						ms = 0;
						for (h = 0; h < 3; h++) {
                  			is[h] = id[h]*(Refine[h]+1) + imins[h]; 
		  					is[h] += (com->facedim == h ? (centered[le][h] == 1 ? 0 : (1-side)*(Refine[h]+1)) : 0);
		  					n += id[h]*stridecr[h];
		  					ms += is[h]*strides[h];
						}
	    /* We know average staggered and centered quantities
		over a lower number of fine zones, which allows us
		to avoid an 'ExecCommSame' communication prior to
		this one. */
						for (iid[2] = 0; iid[2] <= Refine[2]*centered[le][2]; iid[2]++) { 
                  			coef2 =  (1+(1-centered[le][2])*Refine[2]);
		  					for (iid[1] = 0; iid[1] <= Refine[1]*centered[le][1]; iid[1]++) { 
                    			coef1 = coef2 * (1+(1-centered[le][1])*Refine[1]);
	            				for (iid[0] = 0; iid[0] <= Refine[0]*centered[le][0]; iid[0]++) {
                      				coef0 = coef1 * (1+(1-centered[le][0])*Refine[0]);
		      						src = source[le][ms+iid[0]*strides[0]+ iid[1]*strides[1] + iid[2]*strides[2]];
		      						com->buffer[le*com->size+n] += src*coef0/(real)(folded);
		    					}
		  					}
						}
	      			} // i
	    		} // j
	  		} // k 
		} //le
      }
    }
    com = com->next; 
  }
  ExecComm (lev,lev-1,MEAN,nvar,fieldtype);
}

void ExecCommDownFlux (long lev) {
#ifndef GPU
  /* Execute communications from level lev to coarser level lev-1, of
     type FLUX */
  jCommunicator *com;
  long nvar, dim, side, i, imin[3], imax[3], dr[3];
  long size[3], csized[3], csizes[3];
  long dimp1, dimp2, n, le, j, m, buf_size;
  real *source;
  FluidPatch *fluid;
  com = ComListFlux;
  nvar = 1+2*NDIM;		/* Density (or mass flux) + velocities (or momentum flux) */
#ifdef ADIABATIC
  nvar+=1;	/* Energy flux */
#endif
  while (com != NULL) {
    if ((com->dest_level == lev-1) && (com->src_level == lev) && (com->type == FLUX)) {
      if (com->CPU_src == CPU_Rank) {
		dim  = com->facedim;
		side = com->faceside;
		for (i = 0; i < 3; i++) {
		  dr[i] = (Refine[i] ? 2:1);
		  imin[i] = com->imin_src[i];
		  imax[i] = com->imax_src[i];
		  size[i] = com->srcg->ncell[i];
		  csizes[i] = com->imax_src[i]-com->imin_src[i];
		  csized[i] = com->imax_dest[i]-com->imin_dest[i];
		  if ((csizes[i]/dr[i] != csized[i]) && (i != dim))
		    prs_error ("Flux communicator size internal error.");
		}
		fluid = com->srcg->fluid;
		while (fluid != NULL) {
			if(fluid->FluidRank == Current_Fluid){	
		  		source = fluid->Fluxes[dim][side];
			}
		  fluid = fluid->next;
		}
		dimp1 = (dim == 0);
		dimp2 = 2 - (dim == 2); //dim is the dimension orthogonal to the face so we iterate over the other two
		buf_size = csized[dimp1]*csized[dimp2];
		for (i = 0; i < buf_size*nvar; i++)
		  com->buffer[i] = 0.0;
		for (le = 0; le < nvar; le++) {
		  for (i = imin[dimp1]; i < imax[dimp1]; i++) {
		    for (j = imin[dimp2]; j < imax[dimp2]; j++) {
		      m = (i-Nghost[dimp1])+(j-Nghost[dimp2])*size[dimp1];
		      n = (i-imin[dimp1])/dr[dimp1]+(j-imin[dimp2])/dr[dimp2]*csized[dimp1];
		      com->buffer[n+le*buf_size] += source[m+le*size[dimp1]*size[dimp2]];
		    }
		  }
		}
      }
    }
    com = com->next;
  }
  ExecCommFlux (lev);
#else
  ExecCommDownFlux_gpu(lev);
#endif
}

