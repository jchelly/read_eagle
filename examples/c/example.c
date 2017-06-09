#include <stdlib.h>
#include <stdio.h>
#include "read_eagle.h"

int main(int argc, char *argv[])
{
  int i, n;
  float *pos;
  long long *ids;

  abort_on_error(0);

  /* Open the snapshot */
  EagleSnapshot *snap = open_snapshot("./snapshot_002_z009p993/snap_002_z009p993.0.hdf5");  
  if(!snap)
    {
      printf("open_snapshot failed!\n");
      printf("Reason: %s\n", get_error());
      exit(1);
    }

  /* Write out names of available datasets */
  char buf[200];
  int itype, iset, nset;
  for(itype=0;itype<6;itype+=1)
    {
      nset = get_dataset_count(snap, itype);
      for(iset=0;iset<nset;iset+=1)
	{
	  get_dataset_name(snap, itype, iset, buf, 200);
	  printf("# Type %d has dataset: %s\n", itype, buf);
	}
    }

  /* Choose region to read */
  select_region(snap, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0);

  /* Find out how many particles we're going to get back */
  n = count_particles(snap, 0);

  /* Read gas particle positions */
  pos = malloc(sizeof(float)*3*n);
  if(read_dataset_float(snap, 0, "Coordinates", pos, n*3) < 0)
    {
      printf("read_dataset failed!\n");
      printf("Reason: %s\n", get_error());
      exit(1);
    }

  /* Read gas particle IDs */
  ids = malloc(sizeof(long long)*n);
  read_dataset_long_long(snap, 0, "ParticleIDs", ids, n);

  /* Write out the particles */
  for(i=0;i<n;i++)
    printf("%i %14.6f %14.6f %14.6f\n", (int) ids[i], pos[3*i+0], pos[3*i+1], pos[3*i+2]);

  /* Deallocate hash table etc */
  close_snapshot(snap);
}
