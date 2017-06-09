#include <stdio.h>
#include <string.h>
#include "read_eagle.h"

#ifdef CMAKE_F90_NAMES
/* If we're using cmake, this header specifies the F90 name mangling scheme */
#include "FC.h"
#else
/* Otherwise just assume we need to add an underscore to make functions fortran callable */
#define FC_GLOBAL(x,y) x##_
#endif

/*
  Fortran wrappers for read_eagle C functions.

  Intended to be called via read_eagle_fortran.f90.

  Here we assume that strings have been null terminated
  by the fortran caller. Everything else is passed by
  reference.

  Underscores are removed from names to reduce name mangling
  problems.
*/

void FC_GLOBAL(geterrorf, GETERRORF)(int *len, char *str)
{
  strncpy(str, get_error(), (size_t) *len);
}

void FC_GLOBAL(abortonerrorf, ABORTONERRORF)(int *flag)
{
  abort_on_error(*flag);
}

void FC_GLOBAL(opensnapshotf, OPENSNAPSHOTF)(EagleSnapshot **snap, char *fname, double *boxsize,
					     long long *numpart_total, int *numfiles, int *hashbits)
{
  int i;

  *snap = open_snapshot(fname);
  if(*snap)
    {
      *boxsize = (*snap)->boxsize;
      *numfiles = (*snap)->numfiles;
      *hashbits = (*snap)->hashbits;
      for(i=0;i<6;i+=1)
	numpart_total[i] = (*snap)->numpart_total[i];
    }
}

void FC_GLOBAL(closesnapshotf, CLOSESNAPSHOTF)(EagleSnapshot **snap)
{
  close_snapshot(*snap);
}

void FC_GLOBAL(selectregionf, SELECTREGIONF)(EagleSnapshot **snap, 
					     double *xmin, double *xmax,
					     double *ymin, double *ymax,
					     double *zmin, double *zmax)
{
  select_region(*snap,
		*xmin, *xmax,
		*ymin, *ymax,
		*zmin, *zmax);
}

void FC_GLOBAL(selectgridcellsf, SELECTGRIDCELLSF)(EagleSnapshot **snap, 
						   int *ixmin, int *ixmax,
						   int *iymin, int *iymax,
						   int *izmin, int *izmax)
{
  select_grid_cells(*snap,
		    *ixmin, *ixmax,
		    *iymin, *iymax,
		    *izmin, *izmax);
}

void FC_GLOBAL(setsamplingratef, SETSAMPLINGRATEF)(EagleSnapshot **snap, 
		       double *rate)
{
  set_sampling_rate(*snap, *rate);
}

void FC_GLOBAL(clearselectionf, CLEARSELECTIONF)(EagleSnapshot **snap)
{
  clear_selection(*snap);
}

void FC_GLOBAL(countparticlesf, COUNTPARTICLESF)(long long *n, EagleSnapshot **snap, int *itype)
{
  *n = count_particles(*snap, *itype);
}

void FC_GLOBAL(getparticlelocationsf, GETPARTICLELOCATIONSF)(long long *n, EagleSnapshot **snap, int *itype,
							     int *file_index, int *file_offset, long long *nmax)
{
  *n = get_particle_locations(*snap, *itype, file_index, file_offset, 
			      (size_t) *nmax);
}

void FC_GLOBAL(readdatasetf, READDATASETF)(long long *nread, EagleSnapshot **snap, int *itype, int *typecode, void *buf, 
					   long long *n, char *name)
{
  hid_t dtype_id;
  if(*typecode==0)
    dtype_id = H5T_NATIVE_INT;
  else if(*typecode==1)
    dtype_id = H5T_NATIVE_LLONG;
  else if(*typecode==2)
    dtype_id = H5T_NATIVE_FLOAT;
  else if(*typecode==3)
    dtype_id = H5T_NATIVE_DOUBLE;
  
  *nread = read_dataset(*snap, *itype, name, dtype_id, buf, *n);
}

void FC_GLOBAL(splitselectionf, SPLITSELECTIONF)(int *ret, EagleSnapshot **snap, int *ThisTask, int *NTask)
{
  *ret = split_selection(*snap, *ThisTask, *NTask);
}
