#include "fargo3d.h"

//extern long cpugrid_number;

void release_com(CommHash *hash){
  jCommunicator *prev, *next, *com;
  CommHash *hprev, *hnext;
  com = hash->com;
  prev = com->prev;
  next = com->next;
  /* If that communicator is not at the
				   begining of the list : */
  if (prev != NULL)
    prev->next = next;
  /* If that communicator is not at the
				   end of the list : */
  if (next != NULL)
    next->prev = prev;
  hprev = hash->prev;
  hnext = hash->next;
  if (hnext != NULL)
    hnext->prev = hprev;
  if (hprev != NULL)
    hprev->next = hnext;
  /* If the communicator was the first
				   one do not forget to update the
				   global pointer to the chained
				   list */
  if (com == ComListGhost)
    ComListGhost = next;
  /* If the hash entry was the first one
				   update the global pointer in the
				   corresponding array */
  if (hprev == NULL)
    CommHashSrc[hash->com->nb_src] = hnext;
  free(com);
  free(hash);
}

void DestroyHash(){
  CommHash *hash, *next;
  long i;
  for (i = 0; i < cpugrid_number; i++)
  {
    hash = CommHashSrc[i];
    while (hash != NULL)
    {
      next = hash->next;
      free(hash);
      hash = next;
    }
    CommHashSrc[i] = NULL;
  }
}

void DestroyCom(jCommunicator **start){
  jCommunicator *com, *next;
  com = *start;
  while (com != NULL)
  {
    next = com->next;
    free(com);
    *start = com = next;
  }
}

long StripCommunicators()
{
  jCommunicator *com1, *com2;
  CommHash *hash, *next;
  long min1[3], max1[3], min2[3], max2[3], i, stripped = 0;
  boolean included;
  pInfo("Stripping communicators...\n");
  com1 = ComListGhost;
  while (com1 != NULL)
  {
    hash = CommHashSrc[com1->nb_src];
    while (hash != NULL)
    {
      next = hash->next;
      com2 = hash->com;
      if ((com1->nb_dest == com2->nb_dest) &&
          (com1 != com2))
      {
        for (i = 0; i < NDIM; i++)
        {
          min1[i] = com1->Imin[i];
          min2[i] = com2->Imin[i];
          max1[i] = com1->Imax[i];
          max2[i] = com2->Imax[i];
        }
        included = TRUE;
        for (i = 0; i < NDIM; i++)
          if ((min2[i] < min1[i]) || (max2[i] > max1[i]))
            included = NO;
        if (included)
        {
          release_com(hash);
          stripped++;
        }
      }
      hash = next;
    }
    com1 = com1->next;
  }
  pInfo("%ld communicators stripped\n", stripped);
  return stripped;
}
