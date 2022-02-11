#include "fargo3d.h"

extern MPI_Request *Req;

/* Before calling the following function, one needs to fill the
   buffers of the corresponding communicators */
void BuildCommDest(FluidPatch *fluid, int *fieldtype, int nvar)
{
  int i, field, le = 0;
  while (fluid != NULL)
  {
    if (fluid->FluidRank == Current_Fluid)
    {
      for (i = 0; i < nvar; i++)
      {
        field = fieldtype[i];
        centered[le][0] = centered[le][1] = centered[le][2] = 1;
        staggered[le] = 0;
        stdim[le] = 0;
        if (field == _Density_)
        {
          JSINPUT(fluid->Density);
          JSOUTPUT(fluid->Density);
        }
        else if (field == _Energy_)
        {
          JSINPUT(fluid->Energy);
          JSOUTPUT(fluid->Energy);
        }
#ifdef X
        else if (field == _Vx_)
        {
          JVINPUT(fluid->Velocity, _X_);
          JVOUTPUT(fluid->Velocity, _X_);
          centered[le][_X_] = 0;
          staggered[le] = 1;
          stdim[le] = _X_;
        }
#endif
#ifdef Y
        else if (field == _Vy_)
        {
          JVINPUT(fluid->Velocity, _Y_);
          JVOUTPUT(fluid->Velocity, _Y_);
          centered[le][_Y_] = 0;
          staggered[le] = 1;
          stdim[le] = _Y_;
        }
#endif
#ifdef Z
        else if (field == _Vz_)
        {
          JVINPUT(fluid->Velocity, _Z_);
          JVOUTPUT(fluid->Velocity, _Z_);
          centered[le][_Z_] = 0;
          staggered[le] = 1;
          stdim[le] = _Z_;
        }
#endif
        dest[le++] = fluid->Ptr[field];
      }
    }
    fluid = fluid->next;
  }
}

void ExecComm(long levsrc, long levdest, long type, long nvar, int *fieldtype, int facedim)
{
  jCommunicator *com = NULL, *start = NULL;
  long i, nbreq = 0, le, imin[3], imax[3], stride[3];
  long comp[20], j, k, m, n, d, sqz[3], ii[3], h;
  FluidPatch *fluid;
  MPI_Status stat;
  int field, displacement=0;
  tGrid_CPU *desc;
  real s, vp;
  if (type == GHOST)
    start = com = ComListGhost;
  if (type == MEAN)
    start = com = ComListMean;

  while (com != NULL)
  { /* We initiate non-blocking communications */
    if ((com->dest_level == levdest) && (com->src_level == levsrc) && (com->type == type) && (com->facedim == facedim))
    {
      if (com->CPU_src != com->CPU_dest)
      {
        if (com->CPU_src == CPU_Rank)
        {
          MPI_Isend(com->buffer, com->size * nvar, MPI_DOUBLE, com->CPU_dest, com->tag, DomainComm, Req + nbreq);
          nbreq++;
        }
        if (com->CPU_dest == CPU_Rank)
        {
          MPI_Irecv(com->buffer, com->size * nvar, MPI_DOUBLE, com->CPU_src, com->tag, DomainComm, Req + nbreq);
          nbreq++;
        }
      }
    }
    com = com->next;
  }
  for (i = 0; i < nbreq; i++)
    MPI_Wait(Req + i, &stat);
  com = start;

  while (com != NULL)
  { /* We use the buffers to update the
				   ghost or "mean" values */
    if ((com->dest_level == levdest) && (com->src_level == levsrc) && (com->type == type) && (com->facedim == facedim))
    {
      if (com->CPU_dest == CPU_Rank)
      {
        le = 0;
        fluid = com->destg->fluid;
        BuildCommDest(fluid, fieldtype, nvar);
        desc = com->destg;
        for (le = 0; le < 3; le++)
        {
          imin[le] = com->imin_dest[le];
          imax[le] = com->imax_dest[le];
          stride[le] = com->destg->stride[le];
        }
        n = 0;
        for (le = 0; le < nvar; le++)
        {
          for (k = imin[2]; k < imax[2]; k++)
          {
            for (j = imin[1]; j < imax[1]; j++)
            {
              for (i = imin[0]; i < imax[0]; i++)
              {
                ii[0] = i;
                ii[1] = j;
                ii[2] = k;
                /* If the quantity to deploy is a staggered vector
		            perpendicular to the communicator and on a lower
		            face of the destination grid, we move it one cell
		            to the right along the direction of the staggered
		            vector, because this is where the source quantity
		            was considered (there is no position mismatch
		            between source and destination, and the result at a
		            given mesh location should be independent of which
		            kind of communicator (direction X Y or Z) is used
		            to send the data). See matching comment in
		            J_comm_fill.c */
                displacement = 0;
#if COMM != ASYMMETRIC
                if (levdest > levsrc)
                { //GHOST CASE
#if COMM == INTERNAL
                  displacement = (1 - com->faceside) * (1 - centered[le][com->facedim]);
#endif
#if COMM == EXTERNAL
                  displacement = com->faceside * (1 - centered[le][com->facedim]);
#endif
                }
#endif
                if (levdest < levsrc)
                { //MEAN CASE
                  displacement = (1 - com->faceside) * (1 - centered[le][com->facedim]);
                }
                ii[com->facedim] += displacement;
                m = 0;
                for (h = 0; h < 3; h++)
                  m += ii[h] * stride[h];
                s = com->buffer[n++];
#ifdef COMMADAPT
                if (levdest > levsrc && centered[le][_X_] == 0)
                {
                  vp = desc->Ymed[ii[1]] * sin(desc->Zmed[ii[2]]) * OMEGAFRAME;
                  s = s / sqrt(desc->Ymed[ii[1]]);
                  s = s - vp;
                }
#endif
                if (levdest >= levsrc || staggered[le] == 0 ||
                    centered[le][com->facedim] == 0)
                {
                  dest[le][m] = s;
                }
                else if ((ii[stdim[le]] > imin[stdim[le]]) ||
                         (com->srcg->iface[stdim[le]][INF] != -1))
                {
                  dest[le][m] = s;
                }
              }
            }
          }
        }
        if (com->size * nvar != n)
          prs_error("Internal error: communicator size mismatch error");
      }
    }
    com = com->next;
  }
}

void ExecCommFlux(long levsrc)
{
  jCommunicator *com;
  long nbreq = 0, nvar, i;
  long levdest;
  MPI_Status stat;
  levdest = levsrc - 1;
  nvar = 1 + 2 * NDIM; /* Mass + Momentum flux*/
  /* Two momenta flavors per direction (plus and minus) */
#ifdef ADIABATIC
  nvar += 1; /* Energy flux */
#endif
  //printf("Executing Flux comms\n");
  com = ComListFlux;
  while (com != NULL)
  { /* We initiate non-blocking communications */
    if ((com->dest_level == levdest) && (com->src_level == levsrc) && (com->type == FLUX))
    {
      if (com->CPU_src != com->CPU_dest)
      {
        if (com->CPU_src == CPU_Rank)
        {
          MPI_Isend(com->buffer, com->size * nvar, MPI_DOUBLE, com->CPU_dest, com->tag, DomainComm, Req + nbreq);
          nbreq++;
        }
        if (com->CPU_dest == CPU_Rank)
        {
          MPI_Irecv(com->buffer, com->size * nvar, MPI_DOUBLE, com->CPU_src, com->tag, DomainComm, Req + nbreq);
          nbreq++;
        }
      }
    }
    com = com->next;
  }
  for (i = 0; i < nbreq; i++)
    MPI_Wait(Req + i, &stat);
}
