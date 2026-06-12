#define __NOPROTO
#define __GPU

#include "fargo3d.h"
#include "J_jupiter.h"

#define THREADS 128
__device__ real DvanLeer (real slope1, real slope2)
{
  if (slope1 * slope2 <= 0.0)
    return 0;
  return (2.*slope1*slope2/(slope1+slope2));
}

__global__ void kernel_UPLIL (real **source, real *ymin, real *zmin, real OF, struct parms srcParms, real **buffers, struct gpucomm *CommArray, int *centered, int NDIM, int fields){

  int bufferid = blockIdx.x / fields;
  struct gpucomm com = CommArray[bufferid];
  real *buffer = buffers[bufferid];
  int le = blockIdx.x - fields*bufferid;
  int index, pitch, pitchd, stride, strided, size, sized, faceside, facedim;
  int id[3], is[3], icd[3], ics[3], imin[3], dnd[3], dns[3], cornerd[3], corners[3];
  int j, k, h, ms, md, msl;
  real s, src, srcp, srcm, slope,frac,value_i[9],value_j[3],ymed,vp,sr;
  faceside = com.faceside;
  facedim =  com.facedim;
  size = com.size;
  sized = com.sizeD/sizeof(real);
  pitch = com.Pitch;
  stride = com.Stride;
  pitchd = com.pitchD/sizeof(real);
  strided = com.strideD/sizeof(real);
  imin[0] = com.xmin;
  imin[1] = com.ymin;
  imin[2] = com.zmin;
  cornerd[0] = com.cxmin;
  cornerd[1] = com.cymin;
  cornerd[2] = com.czmin;
  dnd[0] = com.dx;
  dnd[1] = com.dy;
  dnd[2] = com.dz;
  corners[0] = srcParms.cxmin;
  corners[1] = srcParms.cymin;
  corners[2] = srcParms.czmin;
  dns[0] = srcParms.dx;
  dns[1] = srcParms.dy;
  dns[2] = srcParms.dz;
 
  
  for(index = threadIdx.x; index < size; index += blockDim.x){
  	id[2] = index/stride;
	id[1] = (index - id[2]*stride)/pitch;
	id[0] = index - id[2]*stride - id[1]*pitch;    
  	for (h = 0; h < 3; h++){
#if COMM == INTERNAL
	  icd[h] = ((id[h]+imin[h])*dnd[h]+(dnd[h]>>1)*centered[le*3+h]+(facedim == h ? (centered[le*3+h] == 1 ? 0 : (1-faceside)*dnd[h]) : 0)+cornerd[h]);// c like 'canvas'
#endif
#if COMM == EXTERNAL
	  icd[h] = ((id[h]+imin[h])*dnd[h]+(dnd[h]>>1)*centered[le*3+h]+(facedim == h ? (centered[le*3+h] == 1 ? 0 : faceside*dnd[h]) : 0)+cornerd[h]);
#endif
#if COMM == ASYMETRIC
	  icd[h] = (id[h]+imin[h])*dnd[h]+(dnd[h]>>1)*centered[le*3+h]+cornerd[h];
#endif
	  is[h] = (icd[h]-corners[h])/dns[h];
	  ics[h] = (is[h]*dns[h]+(dns[h]>>1)*centered[le*3+h]+corners[h]);// c like 'canvas'
	}
	ms = is[0] + is[1]*srcParms.pitch + is[2]*srcParms.stride;
	frac = (double)(icd[0]-ics[0])/(double)(dns[0]);
	for (k=-(NDIM>2); k<=(NDIM>2); k++) {//-1 to +1 if 3D, 0 otherwise.
	  for (j=-(NDIM>1); j<=(NDIM>1); j++) {//-1 to +1 if at least 2D, 0 otherwise.
	    msl = ms+j*srcParms.pitch+k*srcParms.stride; //index on source mesh, local, with offset
	    src = source[le][msl];
	    srcp= source[le][msl+1];
	    srcm= source[le][msl-1];
#ifdef COMMADAPT
	     if(centered[le*3] == 0){
	       ymed = 0.5*(ymin[is[1]+j]+ymin[is[1]+j+1]);
	       vp = ymed*sin(0.5*(zmin[is[2]+k]+zmin[is[2]+k+1]))*OF;
	       sr = sqrt(ymed);
	       src = (src +vp)*sr;
	       srcp= (srcp+vp)*sr;
	       srcm= (srcm+vp)*sr;
	     }
#endif
	     
#ifdef NONMONOTONIC
	     if(frac>=0)
	       value_i[j+1+3*(k+1)]=(1-frac)*src+frac*srcp;
	     else
	       value_i[j+1+3*(k+1)]=(1+frac)*src-frac*srcm;
#else
	     slope = DvanLeer(srcp-src,src-srcm); //Do not divide by mesh size
	     value_i[j+1+3*(k+1)] = src+slope*frac;
#endif
	  }
	}
	s = value_i[4];
	if (NDIM > 1) {
	  frac = (double)(icd[1]-ics[1])/(double)(dns[1]);
	  for (k=-(NDIM>2); k<=(NDIM>2); k++) {//-1 to +1 if 3D, 0 otherwise.
	    src = value_i[1+3*(k+1)];
	    srcp= value_i[2+3*(k+1)];
	    srcm= value_i[  3*(k+1)];
#ifdef NONMONOTONIC
	     if(frac>=0)
	       value_j[k+1]=(1-frac)*src+frac*srcp;
	     else
	       value_j[k+1]=(1+frac)*src-frac*srcm;
#else
	    slope = DvanLeer(srcp-src,src-srcm);
	    value_j[k+1] = src+slope*frac;
#endif
	  }
	  s = value_j[1];
	}
	if (NDIM > 2) {
	  frac = (double)(icd[2]-ics[2])/(double)(dns[2]);
	  src = value_j[1];
	  srcp= value_j[2];
	  srcm= value_j[0];
#ifdef NONMONOTONIC
	     if(frac>=0)
	       s=(1-frac)*src+frac*srcp;
	     else
	       s=(1+frac)*src-frac*srcm;
#else
	    slope = DvanLeer(srcp-src,src-srcm);
	    s = src+slope*frac;
#endif
	}
	md = id[0]+id[1]*pitchd+id[2]*strided;
	buffer[md+sized*le] = s;
	
 }
 
    
}

__global__ void kernel_DOWNMEAN (real **source, struct parms srcParms, real **buffers, struct gpucommDown *CommArray, int *centered, int fields){

  int bufferid = blockIdx.x / fields;
  struct gpucommDown com = CommArray[bufferid];
  real *buffer = buffers[bufferid];
  int le = blockIdx.x - fields*bufferid;
  int index, pitch, pitchd, stride, strided, size, sized, faceside, facedim;
  int id[3], is[3], imin[3], Refine[3], folded;
  int i, j, k, h, ms, md;
  real s, src, coef2, coef1, coef0;
  faceside = com.faceside;
  facedim =  com.facedim;
  size = com.size;
  sized = com.sizeD/sizeof(real);
  pitch = com.Pitch;
  stride = com.Stride;
  pitchd = com.pitchD/sizeof(real);
  strided = com.strideD/sizeof(real);
  imin[0] = com.xmin;
  imin[1] = com.ymin;
  imin[2] = com.zmin;
  Refine[0] = com.xrefine;
  Refine[1] = com.yrefine;
  Refine[2] = com.zrefine;
  
  folded = (Refine[0] ? 2:1) * (Refine[1] ? 2:1) * (Refine[2] ? 2:1);
    
  for(index = threadIdx.x; index < size; index += blockDim.x){
  	id[2] = index/stride;
	id[1] = (index - id[2]*stride)/pitch;
	id[0] = index - id[2]*stride - id[1]*pitch;    
  	for (h = 0; h < 3; h++) {
	  is[h] = id[h]*(Refine[h]+1) + imin[h]; 
	  is[h] += (facedim == h ? (centered[le*3+h] == 1 ? 0 : (1-faceside)*(Refine[h]+1)) : 0);
	}
	ms = is[0] + is[1]*srcParms.pitch + is[2]*srcParms.stride;
	s=0;

	for (k = 0; k <= Refine[2]*centered[le*3+2]; k++) { 
            coef2 =  (1+(1-centered[le*3+2])*Refine[2]);
	    for (j = 0; j <= Refine[1]*centered[le*3+1]; j++) { 
                coef1 = coef2 * (1+(1-centered[le*3+1])*Refine[1]);
	        for (i = 0; i <= Refine[0]*centered[le*3]; i++) {
                    coef0 = coef1 * (1+(1-centered[le*3])*Refine[0]);
	            src = source[le][ms+i+j*srcParms.pitch + k*srcParms.stride];
		    s += src*coef0/(real)(folded);
		}
	     }
	}
	md = id[0]+id[1]*pitchd+id[2]*strided;
	buffer[md+sized*le] = s;
 }
 
    
}


__global__ void kernel_DOWNFLUX (real **source, struct parms srcParms, real **buffers, struct gpucommFlux *CommArray, int fields, int cpu){

  int bufferid = blockIdx.x / fields;
  struct gpucommFlux com = CommArray[bufferid];
  real *buffer = buffers[bufferid];
  int le = blockIdx.x - fields*bufferid;
  int index, pitchd, strided, size, sized, side, dim, dim1, dim2;
  int id[3], is[3], Refine[3], sizes[3], imin[3], n1;
  int i, j, h, ms, md;
  real s=0;
  side = com.faceside;
  dim =  com.facedim;
  dim1 = (dim == 0);
  dim2 = 2 - (dim == 2);
  size = com.size;
  sized = com.sizeD/sizeof(real);
  n1 = com.n1;
  pitchd = com.pitchD/sizeof(real);
  strided = com.strideD/sizeof(real);
  Refine[0] = com.xrefine+1;
  Refine[1] = com.yrefine+1;
  Refine[2] = com.zrefine+1;
  sizes[0] = srcParms.sizex;
  sizes[1] = srcParms.sizey;
  sizes[2] = srcParms.sizez;
  imin[0] = com.xmin;
  imin[1] = com.ymin;
  imin[2] = com.zmin;
    
  for(index = threadIdx.x; index < size; index += blockDim.x){
    s=0;
    id[dim2] = index/n1;
    id[dim1] = index - id[dim2]*n1;
    id[dim]= 0;
    for(h=0;h<3;h++)
      is[h] = (id[h]*Refine[h])+imin[h];
    for(j=0; j < Refine[dim2]; j++){
      for(i=0; i < Refine[dim1]; i++){
	ms = (is[dim1]+i) + (is[dim2]+j) * sizes[dim1];
	s += source[dim*2+side][ms+le*sizes[dim1]*sizes[dim2]];
      }
    }
    md = id[0]+id[1]*pitchd+id[2]*strided;
    buffer[md+sized*le] = s;
  }
       
}

extern "C"
void UPLIL (tGrid_CPU *grid, int fields)
{
 if(grid->src.nbCommsUp !=0){
   kernel_UPLIL<<<grid->src.nbCommsUp*fields,THREADS>>>(grid->source, grid->Ymin_d, grid->Zmin_d, OMEGAFRAME, grid->gpuparms, grid->src.BuffersUp, grid->src.ParmsUp, grid->centered_gpu, NDIM, fields);
   cudaThreadSynchronize();
 }
}		

extern "C"
void DOWNMEAN (tGrid_CPU *grid, int fields)
{
  if(grid->src.nbCommsDown !=0){
    kernel_DOWNMEAN<<<grid->src.nbCommsDown*fields,THREADS>>>(grid->source,  grid->gpuparms, grid->src.BuffersDown, grid->src.ParmsDown, grid->centered_gpu, fields);
    cudaThreadSynchronize();
  }
}		

extern "C"
void DOWNFLUX (tGrid_CPU *grid, int fields)
{
  if(grid->src.nbCommsFlux !=0){
    kernel_DOWNFLUX<<<grid->src.nbCommsFlux*fields,THREADS>>>(grid->Fluid->FluxesGPU, grid->gpuparms, grid->src.BuffersFlux, grid->src.ParmsFlux, fields, CPU_Rank);
    cudaThreadSynchronize();
  }
}
