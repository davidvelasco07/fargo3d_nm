#include "fargo3d.h"
#include "J_jupiter.h"

int Inside_sphere(int i,int j,int k,real xb,real yb,real zb){
  real sm2=SMOOTHING*SMOOTHING;
  real r2,x,y,z;
  x=XC-xb;
  y=YC-yb;
  z=ZC-zb;
  r2=x*x+y*y+z*z;
  return (r2<=sm2);
}

void Count_Luminous_cells(){
  int i,j,k,ii,jj,kk;
  real xb,yb,zb;
  tGrid_CPU *item;
  item = Current_Jupiter_Patch;
  for(kk=-1;kk<=1;kk++){
    for(jj=-1;jj<=1;jj++){
      for(ii=-1;ii<=1;ii++){
	xb=(abs(ii)>0?2:0)
	  *(ii<0?item->Edges[0][NGHX]:item->Edges[0][Nx+NGHX+1])
	  +XBODY*(ii==0?1:-1);
	yb=(abs(jj)>0?2:0)
	  *(jj<0?item->Edges[1][NGHY]:item->Edges[1][Ny+NGHY+1])
	  +YBODY*(jj==0?1:-1);
	zb=(abs(kk)>0?2:0)
	  *(kk<0?item->Edges[2][NGHZ]:item->Edges[2][Nz+NGHZ+1])
	  +ZBODY*(kk==0?1:-1);
	for(k=NGHZ;k<Nz+NGHZ;k++)
	  for(j=NGHY;j<Ny+NGHY;j++)
	    for(i=NGHX;i<NX+NGHX;i++)
	      item->LumCells+=Inside_sphere(i,j,k,xb,yb,zb);
	
      }
    }
  }
}


