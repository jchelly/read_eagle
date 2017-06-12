
Read routines for Peano-Hilbert key sorted Eagle snapshots
----------------------------------------------------------

These can be used if Gadget was run with the -DEAGLE_SORT_OUTPUT
flag enabled. They provide a quick way to read in spatial regions
without having to read all the data or even open all of the files.

It works by splitting the simulation volume into a grid and ensuring
that particles in the same grid cell are stored consecutively in the 
snapshot files. Some extra datasets are added to the snapshot which
specify which grid cells are stored in which files and the location
of the start of each cell in the file.

The procedure to read in data is as follows:

- Open the snapshot by specifying the name of one file
- Flag the grid cells to be read in by calling select_region() one or
  more times.
- Call read_dataset once for each quantity to be read, specifying
  which particle type to read and the HDF5 dataset name

You can then either close the snapshot or call clear_selection() to
read in a different region.

Any of the datasets in the PartType groups can be read in. The code
makes no assumptions about what datasets are present in the snapshot.
Reduced precision floating point datasets are read in as 32 bit floats.

The code is written in C but there are interfaces for Fortran
and Python. There's also a pure IDL version since calling the C 
routines from IDL is impractical due to library version issues.


Installing the python module
----------------------------

To install the python module:

 - Edit setup.py to specify the location of your HDF5 installation
 - To install to your home directory, run

     python ./setup.py install --user

   or to install elsewhere

     python ./setup.py install --prefix=/path/to/install/dir/

If you use the --prefix option you'll probably need to set the environment
variable PYTHONPATH so that python can find the module.


Compiling the library and examples for C and Fortran
----------------------------------------------------

The library is compiled using cmake. If your HDF5 installation is in
an unusual location you can specify it with cmake's CMAKE_PREFIX_PATH
option. E.g.:

  mkdir build
  cd build
  cmake .. -DCMAKE_PREFIX_PATH=/path/to/hdf5/installation
  make

The path you pass in should contain the HDF5 lib and include directories.
This will produce the following files:

  build/lib/libread_eagle.so   - the shared library for reading Eagle snapshots
  build/include/read_eagle.h   - header file needed to call the library from C
  build/include/read_eagle.mod - module needed to call the library from Fortran
  build/bin/*                  - compiled example programs in C and Fortran


Example usage in C:
-------------------

  #include "read_eagle.h"

  EagleSnapshot *snap;
  float *pos;
  long long *ids;
  int n;

  /* Open snapshot by supplying name of one file */
  snap = open_snapshot("./snap_020.0.hdf5");

  /* Choose region to read by specifying ranges in x, y, z */
  select_region(snap, 4.0, 5.0, 2.0, 3.0, 3.0, 4.0);

  /* Allocate memory */
  n = count_particles(snap, 0);
  pos = malloc(sizeof(float)*3*n);
  ids = malloc(sizeof(long long)*n);

  /* Read one or more datasets */
  n = read_dataset_float(snap, 0, "Coordinates", pos, n*3)
  n = read_dataset_long_long(snap, 0, "ParticleIDs", ids, n)

  /* Close snapshot */
  close_snapshot(snap);


Example usage in Fortran:
-------------------------

  use read_eagle

  type(EagleSnapshot) :: snap
  real, dimension(:,:), allocatable :: pos
  integer*8, dimension(:,:), allocatable :: ids
  integer :: ierr
  integer :: n

  ! Open snapshot by supplying name of one file
  snap = open_snapshot("./snap_020.0.hdf5")

  ! Choose region to read by specifying ranges in x, y, z
  call select_region(snap, 4.0, 5.0, 2.0, 3.0, 3.0, 4.0)

  ! Allocate memory
  n = count_particles(snap, 0);
  allocate(pos(3,n))
  allocate(ids(n))

  ! Read one or more datasets
  n = read_dataset(snap, 0, "Coordinates", pos)
  n = read_dataset(snap, 0, "ParticleIDs", ids)

  ! Close snapshot
  call close_snapshot(snap)


Example usage in Python:
------------------------

import read_eagle

snap = read_eagle.EagleSnapshot("./snap_020.0.hdf5")
snap.select_region(4.0, 5.0, 2.0, 3.0, 3.0, 4.0)
pos = snap.read_dataset(0, "Coordinates")
ids = snap.read_dataset(0, "ParticleIDs")
del snap


Example usage in IDL:
---------------------

.run read_eagle.pro
snap = open_snapshot("./snap_020.0.hdf5")
select_region, snap, 4.0, 5.0, 2.0, 3.0, 3.0, 4.0
n = count_particles(snap, 0)
if n gt 0 then begin
  pos = read_dataset(snap, 0, "Coordinates")
  ids = read_dataset(snap, 0, "ParticleIDs")
endif
close_snapshot, snap


Error handling
--------------

In C:

Default behaviour is to abort if an error occurs, e.g. if a file can't
be read or memory allocation fails. If you call abort_on_error(0)
the routines will instead indicate errors with special return values.

open_snapshot   - returns a null pointer on failure
count_particles - returns a negative number on failure
read_dataset    - returns a negative number on failure

In Fortran:

The routines which can fail (open_snapshot, count_particles and
read_dataset) take an optional integer parameter 'iostat'. If this
parameter is not specified the routines will abort on errors. If it
is specified it will return zero on success, non-zero otherwise.

In Python:

If a routine fails it will raise an exception.

In IDL:

The IDL routines always abort in case of errors.



Description of the routines and parameters
------------------------------------------


* open_snapshot
---------------

C      :  EagleSnapshot *snap = open_snapshot(char *fname)
F90    :  snap = open_snapshot(fname [, iostat=...])
Python :  snap = read_eagle.EagleSnapshot(fname)
IDL    :  snap = open_snapshot(fname)

This opens the snapshot which contains the specified file. The returned
value 'snap' must be deallocated using a call to close_snapshot.
Failure to do this in C, Fortran or IDL will cause memory leaks. In 
python the snapshot will be automatically deallocated once there are no
more references to it.

Parameters

  - fname: name of any one file in the snapshot

Return value

  C     : pointer to a newly allocated EagleSnapshot, or NULL on failure
  F90   : an instance of the eaglesnapshot derived type
  Python: an instance of the EagleSnapshot class
  IDL   : a struct containing data necessary to read the snapshot


* close_snapshot
----------------

C       : close_snapshot(EagleSnapshot *snap)
F90     : call close_snapshot(snap)
python  : del snap
IDL     : close_snapshot, snap

Deallocates memory associated with the snap object from open_snapshot.
Not necessary in python - do 'del snap' or just let the snap variable go
out of scope.

Parameters

  - snap: the object returned by open_snapshot


* select_region
---------------

C       : select_region(EagleSnapshot *snap, double xmin, double xmax,
                                             double ymin, double ymax,
                                             double zmin, double zmax)
F90     : call select_region(snap, xmin, xmax, ymin, ymax, zmin, zmax,
                             [, iostat=...])
python  : snap.select_region(xmin, xmax, ymin, ymax, zmin, zmax)
IDL     : select_region, snap, xmin, xmax, ymin, ymax, zmin, zmax

All grid cells overlapping the specified region are flagged to be read
in by subsequent read_dataset calls. You can call select_region multiple
times to make oddly shaped or disjoint selections.

If selected regions overlap or the same region is selected multiple times
particles in these regions will still only be read in once.

Parameters

  - snap: the object returned by open_snapshot
  - xmin: the minimum x coordinate of the region to read
  - xmax: the maximum x coordinate of the region to read
  - ymin: the minimum y coordinate of the region to read
  - ymax: the maximum y coordinate of the region to read
  - zmin: the minimum z coordinate of the region to read
  - zmax: the maximum z coordinate of the region to read

  The coordinates are doubles in C and reals in Fortran.


* select_grid_cells
-------------------

C       : select_grid_cells(EagleSnapshot *snap, int ixmin, int ixmax,
                           		         int iymin, int iymax,
                                                 int izmin, int izmax)
F90     : call select_grid_cells(snap, ixmin, ixmax, iymin, iymax, izmin, izmax,
                                 [, iostat=...])
python  : snap.select_grid_cells(ixmin, ixmax, iymin, iymax, izmin, izmax)
IDL     : Not implemented

All grid cells in the specified range of grid coordinates are flagged to
be read in by subsequent read_dataset calls. You can call select_grid_cells
multiple times to make oddly shaped or disjoint selections.

The coordinates ixmin, ixmax etc are integer coordinates in the hash grid,
starting from zero. The maximum coordinate is (2**hashbits)-1. The value of
hashbits is stored in the EagleSnapshot variable returned by open_snapshot call.

If selected regions overlap or the same region is selected multiple times
particles in these regions will still only be read in once.

Parameters

  - snap: the object returned by open_snapshot
  - ixmin: the minimum x coordinate of the region to read
  - ixmax: the maximum x coordinate of the region to read
  - iymin: the minimum y coordinate of the region to read
  - iymax: the maximum y coordinate of the region to read
  - izmin: the minimum z coordinate of the region to read
  - izmax: the maximum z coordinate of the region to read

  The coordinates are ints in C and default integers in Fortran.


* count_particles
-----------------

C      : int n = count_particles(EagleSnapshot *snap, int itype)
F90    : n = count_particles(snap, itype [, iostat=...])
Python : n = snap.count_particles(itype)
IDL    : n = count_particles(itype)

This returns the number of particles of the specified type which will be
read by the next read_dataset call. Note that only whole grid cells can
be read so some particles outside the selected region may be read in.
These are included in the count.

In C and Fortran this can be used to determine how much memory to allocate
before reading datasets. In python its not usually necessary.

In IDL you need to call this before read_dataset to make sure there are
particles in the selected region - IDL doesn't allow zero size arrays so
read_dataset will abort if there are no particles.


Parameters

  - snap : the object returned by open_snapshot
  - itype: which particle type to count (integer, 0-5)

Return value

  The number of particles to be read in


* get_particle_locations
------------------------

C      : int n = get_particle_locations(EagleSnapshot *snap, int itype,
                                        int *file_index, int *file_offset,
                                        size_t nmax);
F90    : n = get_particle_locations(snap, itype, file_index, file_offset [, iostat=...])
Python : file_index, file_offset = snap.get_particle_locations(itype)
IDL    : Not implemented

This returns two arrays which each have one element for each selected
particle. file_index contains the index of the file each particle is in.
file_offset contains the position in the file, numbering from zero.

In C and Fortran the output arrays must be allocated in advance. Call
count_particles() to determine the necessary size. In C, nmax must be
the size of the file_offset/file_index arrays. This is used to check
the arrays are big enough.

Parameters

  - snap : the object returned by open_snapshot
  - itype: which particle type to count (integer, 0-5)

Return value

  file_index  - integer array with index of the file containing each particle
  file_offset - integer array with position of each particle in its file


* read_dataset
--------------

C       : int read_dataset_int(EagleSnapshot *snap, int itype, char *name, int *buf, size_t n);
          int read_dataset_long_long(EagleSnapshot *snap, int itype, char *name, long long *buf, size_t n);
          int read_dataset_float(EagleSnapshot *snap, int itype, char *name, float *buf, size_t n);
          int read_dataset_double(EagleSnapshot *snap, int itype, char *name, double *buf, size_t n);
Fortran : call read_dataset(snap, itype, name, buf [, iostat=...])
Python  : buf = snap.read_dataset(itype, name)
IDL     : buf = read_dataset(snap, itype, name)

This reads in the specified dataset for all particles of type itype in the
selected region(s). Use repeated calls to read in multiple datasets.

In C the result array 'buf' can be int, long long, float or double. If the
array is not of the same type as the dataset in the file the values will
be converted, which may involve some truncation or rounding.

In Fortran the allowed types for 'buf' are integer*4, integer*8, real*4,
and real*8. Scalar quantities must be read into a 1D array and vectors
must be read into a 2D 3*n array. Type conversion is done in the same way
as in C.

In IDL and Python the type of array you get back reflects the type of the
dataset in the file. The IDL version aborts if there are no particles because
IDL provides no way to return a zero size array section.

Parameters

  - snap : the object returned by open_snapshot
  - itype: which particle type to read (integer, 0-5)
  - name : the HDF5 name of the dataset, relative to the PartType group
  - buf  : array in which to store the result. Must be pre-allocated in
           C or Fortran.
  - n    : size of the array buf (C only, used to check bounds)

Return value

  The number of particles read in


* clear_selection
-----------------

C        : void clear_selection(EagleSnapshot *snap)
Fortran  : call clear_selection(snap)
Python   : snap.clear_selection()
IDL      : clear_selection, snap

Parameters:

  - snap : the object returned by open_snapshot

Return value

  None

This clears the flags which specify which grid cells should be
read in on the next read_dataset() call. If you've already read
in a region and you want to read a different region you should call 
clear_selection() before calling select_region() again.


* split_selection
-----------------

C        : int split_selection(EagleSnapshot *snap, int ThisTask, int NTask)
Fortran  : call split_selection(snap, ThisTask, NTask [, iostat=...])
Python   : snap.split_selection(ThisTask, NTask)
IDL      : Not implemented

Parameters:

  - snap     : the object returned by open_snapshot
  - ThisTask : rank of this process in an MPI program
  - NTask    : number of processes in an MPI program

Return value

  In C or Fortran 0 indicates success, non-zero is failure

This is for use in MPI programs to allow parallel I/O.
When running under MPI all read_eagle routines are collective -
they must be called on all processes.

The procedure to read a region in parallel is as follows:

1. All processes open the snapshot
2. All processes select the SAME region
3. Call split_selection. This causes each processor to select
   a subset of the originally selected particles.
4. Call read_dataset() on all processors to read the required 
   quantities

This results in all of the particles in the specified region 
being read in exactly once with the data spread across the MPI 
processes.
