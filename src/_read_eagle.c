#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NPY_NO_DEPRECATED_API 7
#include "numpy/arrayobject.h"
#include "read_eagle.h"


#if ((PY_MAJOR_VERSION == 2) && (PY_MINOR_VERSION < 7)) || (PY_MAJOR_VERSION < 2)
/* For python < 2.7, we use the CObject API */
#define PTR_FROM_OBJECT(x) PyCObject_AsVoidPtr((x))
#define OBJECT_FROM_PTR(x) PyCObject_FromVoidPtr((x), NULL)
#else
/* For python >= 2.7, we use the Capsule API */
#define OBJECT_FROM_PTR(x) PyCapsule_New((x), NULL, NULL)
#define PTR_FROM_OBJECT(x) PyCapsule_GetPointer((x), NULL)
#endif

/*
  Python interface for read_eagle.c

  Not to be called directly - use the read_eagle.py python
  module.
*/

#define PY_ARRAY_UNIQUE_SYMBOL _READ_EAGLE

/* New python exception for error reporting */
static PyObject *_read_eagleError;

/* Function prototypes */
static PyObject *_read_eagle_open_snapshot(PyObject *self, PyObject *args);
static PyObject *_read_eagle_close_snapshot(PyObject *self, PyObject *args);
static PyObject *_read_eagle_select_region(PyObject *self, PyObject *args);
static PyObject *_read_eagle_select_rotated_region(PyObject *self, PyObject *args);
static PyObject *_read_eagle_select_grid_cells(PyObject *self, PyObject *args);
static PyObject *_read_eagle_set_sampling_rate(PyObject *self, PyObject *args);
static PyObject *_read_eagle_clear_selection(PyObject *self, PyObject *args);
static PyObject *_read_eagle_count_particles(PyObject *self, PyObject *args);
static PyObject *_read_eagle_get_particle_locations(PyObject *self, PyObject *args);
static PyObject *_read_eagle_read_dataset(PyObject *self, PyObject *args);
static PyObject *_read_eagle_read_extra_dataset(PyObject *self, PyObject *args);
static PyObject *_read_eagle_get_dataset_count(PyObject *self, PyObject *args);
static PyObject *_read_eagle_get_dataset_name(PyObject *self, PyObject *args);
static PyObject *_read_eagle_split_selection(PyObject *self, PyObject *args);

/* Method table for the module */
static PyMethodDef _read_eagleMethods[] = {
  {"open_snapshot",   _read_eagle_open_snapshot,   METH_VARARGS, "Open an Eagle snapshot"},
  {"close_snapshot",  _read_eagle_close_snapshot,  METH_VARARGS, "Close an Eagle snapshot"},
  {"select_region",   _read_eagle_select_region,   METH_VARARGS, "Select a region to read in"},
  {"select_rotated_region",   _read_eagle_select_rotated_region,   METH_VARARGS, "Select a non axis aligned region to read in"},
  {"select_grid_cells", _read_eagle_select_grid_cells, METH_VARARGS, "Select grid cells to read in"},
  {"set_sampling_rate", _read_eagle_set_sampling_rate, METH_VARARGS, "Set the sampling rate"},
  {"clear_selection", _read_eagle_clear_selection, METH_VARARGS, "Clear any existing selection"},
  {"count_particles", _read_eagle_count_particles, METH_VARARGS, "Return the number of particles to read"},
  {"get_particle_locations", _read_eagle_get_particle_locations, METH_VARARGS, "Return the locations of particles to read"},
  {"read_dataset",    _read_eagle_read_dataset,    METH_VARARGS, "Read in the selected elements"},
  {"read_extra_dataset",    _read_eagle_read_extra_dataset,    METH_VARARGS, "Read in the selected elements from extra file"},
  {"get_dataset_count", _read_eagle_get_dataset_count, METH_VARARGS, "Return the number of datasets for a type"},
  {"get_dataset_name",  _read_eagle_get_dataset_name,  METH_VARARGS, "Return the name of a dataset"},
  {"split_selection",    _read_eagle_split_selection,    METH_VARARGS, "Split the selection between processors"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

/* Module definition struct for Python 3 only */
#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_read_eagle",
  NULL,
  -1,
  _read_eagleMethods
};
#endif

/* Initialise the module */
#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit__read_eagle(void)
#else
PyMODINIT_FUNC init_read_eagle(void)
#endif
{
  PyObject *m;

  /* Don't abort on error because this would crash the interpreter! */
  abort_on_error(0);

  /* Make sure NumPy is imported and initialised */
  import_array();

#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
  if (m == NULL)
    return NULL;
#else
  m = Py_InitModule("_read_eagle", _read_eagleMethods);
  if (m == NULL)
    return;
#endif

  _read_eagleError = PyErr_NewException("_read_eagle.error", NULL, NULL);
  Py_INCREF(_read_eagleError);
  PyModule_AddObject(m, "error", _read_eagleError);
#if PY_MAJOR_VERSION >= 3
  return m;
#endif
}

/* Open a snapshot */
static PyObject *_read_eagle_open_snapshot(PyObject *self, PyObject *args)
{
  char *fname;
  EagleSnapshot *snap;

  /* Get filename from arguments */
  if (!PyArg_ParseTuple(args, "s", &fname)) return NULL;

  /* Open the file */
  snap = open_snapshot(fname);
  
  /* Raise an exception if anything went wrong */
  if(!snap)
    {
      PyErr_SetString(_read_eagleError, get_error());
      return NULL;
    }

  /* Return the pointer and some info about the snapshot*/
  return Py_BuildValue("NLLLLLLdii", 
		       OBJECT_FROM_PTR(snap),
		       snap->numpart_total[0],
		       snap->numpart_total[1],
		       snap->numpart_total[2],
		       snap->numpart_total[3],
		       snap->numpart_total[4],
		       snap->numpart_total[5],
		       snap->boxsize,
		       snap->numfiles,
		       snap->hashbits);
}

/* Close a snapshot */
static PyObject *_read_eagle_close_snapshot(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;

  /* Get snapshot object from arguments */
  if (!PyArg_ParseTuple(args, "O", &pysnap)) return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);

  /* Close the snapshot */
  close_snapshot(snap);

  Py_INCREF(Py_None);
  return Py_None;
}

/* Choose region to read */
static PyObject *_read_eagle_select_region(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;
  double xmin, xmax;
  double ymin, ymax;
  double zmin, zmax;

  /* Get snapshot and coordinate range from arguments */
  if (!PyArg_ParseTuple(args, "Odddddd", &pysnap,
			&xmin, &xmax, &ymin, &ymax, &zmin, &zmax))
    return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);

  /* Select the cells */
  select_region(snap, xmin, xmax, ymin, ymax, zmin, zmax);

  Py_INCREF(Py_None);
  return Py_None;
}

/* Choose non axis aligned region to read */
static PyObject *_read_eagle_select_rotated_region(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;
  double cx,cy,cz,xx,xy,xz,yx,yy,yz,zx,zy,zz,lx,ly,lz;

  /* Get snapshot and coordinate range from arguments */
  if (!PyArg_ParseTuple(args, "Oddddddddddddddd", &pysnap,
			&cx,&cy,&cz,
			&xx,&xy,&xz,&yx,&yy,&yz,&zx,&zy,&zz,
			&lx,&ly,&lz))
    return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);

  double centre[3] = {cx, cy, cz};
  double xvec[3]   = {xx, xy, xz};
  double yvec[3]   = {yx, yy, yz};
  double zvec[3]   = {zx, zy, zz};
  double length[3] = {lx, ly, lz};

  /* Select the cells */
  select_rotated_region(snap, centre, xvec, yvec, zvec, length);

  Py_INCREF(Py_None);
  return Py_None;
}

/* Choose region to read */
static PyObject *_read_eagle_select_grid_cells(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;
  int ixmin, ixmax;
  int iymin, iymax;
  int izmin, izmax;

  /* Get snapshot and coordinate range from arguments */
  if (!PyArg_ParseTuple(args, "Oiiiiii", &pysnap,
			&ixmin, &ixmax, &iymin, &iymax, &izmin, &izmax))
    return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);

  /* Select the cells */
  select_grid_cells(snap, ixmin, ixmax, iymin, iymax, izmin, izmax);

  Py_INCREF(Py_None);
  return Py_None;
}

/* Set sample rate */
static PyObject *_read_eagle_set_sampling_rate(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;
  double rate;

  /* Get snapshot and sample rate from arguments */
  if (!PyArg_ParseTuple(args, "Od", &pysnap, &rate))
    return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);

  /* Set the sampling rate */
  set_sampling_rate(snap, rate);

  Py_INCREF(Py_None);
  return Py_None;
}


/* Select no region */
static PyObject *_read_eagle_clear_selection(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;

  /* Get snapshot object from arguments */
  if (!PyArg_ParseTuple(args, "O", &pysnap)) return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);
  
  clear_selection(snap);
  
  Py_INCREF(Py_None);
  return Py_None;
}

/* Return the number of particles to read */
static PyObject *_read_eagle_count_particles(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;
  long long n;
  int itype;

  /* Get snapshot object from arguments */
  if (!PyArg_ParseTuple(args, "Oi", &pysnap, &itype)) return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);
  
  /* Get the particle count. n<0 indicates something bad */
  n = count_particles(snap, itype);
  if(n < 0)
    {
      PyErr_SetString(_read_eagleError, get_error());
      return NULL;
    }
  return Py_BuildValue("L", n);
}


/* Return the locations of particles to read */
static PyObject *_read_eagle_get_particle_locations(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;
  long long n;
  int itype;
  int rank;
  npy_intp dims[1];
  PyArrayObject *file_index_array;
  PyArrayObject *file_offset_array;

  /* Get snapshot object from arguments */
  if (!PyArg_ParseTuple(args, "Oi", &pysnap, &itype)) return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);
  
  /* Get the particle count. n<0 indicates something bad */
  n = count_particles(snap, itype);
  if(n < 0)
    {
      PyErr_SetString(_read_eagleError, get_error());
      return NULL;
    }

  /* Construct output arrays */
  rank    = 1;
  dims[0] = n;
  file_index_array  = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_INT);
  file_offset_array = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_INT);
  
  /* Call function to get locations */
  get_particle_locations(snap, itype, 
			 (int *) PyArray_DATA(file_index_array),
			 (int *) PyArray_DATA(file_offset_array),
			 (size_t) n);

  /* Return results */
  return Py_BuildValue("OO", file_index_array, file_offset_array);
}


/* 
   Read in the data

   This needs to behave slightly differently from the C routine.
   We determine the size and type of the result, allocate a 
   new numpy array, read the data, and return the new array.
 */
static PyObject *_read_eagle_read_dataset(PyObject *self, PyObject *args)
{ 
  PyObject      *pysnap;
  EagleSnapshot *snap;
  int itype; 
  char *dset_name;
  long long n;
  TypeCode tc;
  int rank;
  npy_intp dims[2];
  PyArrayObject *result;
  void *data;
  size_t data_size;
  hid_t hdf5_type;

  /* Get parameters from inout arguments */
  if (!PyArg_ParseTuple(args, "Ois", &pysnap, &itype, &dset_name)) return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);
  
  /* Determine number of particles to read */
  n = count_particles(snap, itype);
  if(n < 0)
    {
      PyErr_SetString(_read_eagleError, get_error());
      return NULL;
    }
  
  /* Get type and rank of the dataset and allocate output array */
  if(get_dataset_info(snap, itype, dset_name, &tc, &rank))
    {
      PyErr_SetString(_read_eagleError, get_error());
      return NULL;
    }
  dims[0] = n;
  dims[1] = 3;
  if(tc==t_int)
    {
      hdf5_type = H5T_NATIVE_INT;
      result = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_INT);
    }
  else if(tc==t_long_long)
    {
      hdf5_type = H5T_NATIVE_LLONG;
      result = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_LONGLONG);
    }
  else if(tc==t_float)
    {
      hdf5_type = H5T_NATIVE_FLOAT;
      result = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_FLOAT);
    }
  else
    {
      hdf5_type = H5T_NATIVE_DOUBLE;
      result = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_DOUBLE);
    }

  /* Get a pointer to the data buffer */
  data = PyArray_DATA(result);
  
  /* Read in the data */
  if(rank==1)
    data_size = n;
  else
    data_size = 3*n;
  n = read_dataset(snap, itype, dset_name, hdf5_type, data, data_size);
  if(n < 0)
    {
      /* Read failed so deallocate output array and raise an exception */
      Py_DECREF((PyObject *) result);
      PyErr_SetString(_read_eagleError, get_error());
      return NULL;
    }

  /* Return the result */
  return PyArray_Return(result);
}


/* 
   Read in the data

   This needs to behave slightly differently from the C routine.
   We determine the size and type of the result, allocate a 
   new numpy array, read the data, and return the new array.

   This version is for reading auxiliary files stored in the same
   order as the snapshots.
 */
static PyObject *_read_eagle_read_extra_dataset(PyObject *self, PyObject *args)
{ 
  PyObject      *pysnap;
  EagleSnapshot *snap;
  int itype; 
  char *dset_name;
  long long n;
  TypeCode tc;
  int rank;
  npy_intp dims[2];
  PyArrayObject *result;
  void *data;
  size_t data_size;
  hid_t hdf5_type;
  char *extra_basename;

  /* Get parameters from inout arguments */
  if (!PyArg_ParseTuple(args, "Oiss", &pysnap, &itype, &dset_name, &extra_basename)) return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);
  
  /* Determine number of particles to read */
  n = count_particles(snap, itype);
  if(n < 0)
    {
      PyErr_SetString(_read_eagleError, get_error());
      return NULL;
    }
  
  /* Get type and rank of the dataset and allocate output array */
  if(get_extra_dataset_info(snap, itype, dset_name, &tc, &rank, extra_basename))
    {
      PyErr_SetString(_read_eagleError, get_error());
      return NULL;
    }
  dims[0] = n;
  dims[1] = 3;
  if(tc==t_int)
    {
      hdf5_type = H5T_NATIVE_INT;
      result = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_INT);
    }
  else if(tc==t_long_long)
    {
      hdf5_type = H5T_NATIVE_LLONG;
      result = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_LONGLONG);
    }
  else if(tc==t_float)
    {
      hdf5_type = H5T_NATIVE_FLOAT;
      result = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_FLOAT);
    }
  else
    {
      hdf5_type = H5T_NATIVE_DOUBLE;
      result = (PyArrayObject *) PyArray_SimpleNew(rank, dims, NPY_DOUBLE);
    }

  /* Get a pointer to the data buffer */
  data = PyArray_DATA(result);
  
  /* Read in the data */
  if(rank==1)
    data_size = n;
  else
    data_size = 3*n;
  n = read_extra_dataset(snap, itype, dset_name, hdf5_type, data, data_size, extra_basename);
  if(n < 0)
    {
      /* Read failed so deallocate output array and raise an exception */
      Py_DECREF((PyObject *) result);
      PyErr_SetString(_read_eagleError, get_error());
      return NULL;
    }

  /* Return the result */
  return PyArray_Return(result);
}




/* Return the number of datasets for a particle type */
static PyObject *_read_eagle_get_dataset_count(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;
  int itype;

  /* Get snapshot object and particle type from arguments */
  if (!PyArg_ParseTuple(args, "Oi", &pysnap, &itype)) return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);
  
  /* Range check on itype */
  if(itype < 0 || itype > 5)
    {
      PyErr_SetString(_read_eagleError, "Particle type out of range");
      return NULL;
    }

  /* Return the number of datasets */
  return Py_BuildValue("i", snap->num_datasets[itype]);
}


/* Return the name of a dataset */
static PyObject *_read_eagle_get_dataset_name(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;
  int itype, iset;

  /* Get snapshot object and particle type from arguments */
  if (!PyArg_ParseTuple(args, "Oii", &pysnap, &itype, &iset)) return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);
  
  /* Range check on itype */
  if(itype < 0 || itype > 5)
    {
      PyErr_SetString(_read_eagleError, "Particle type out of range");
      return NULL;
    }

  /* Range check on iset */
  if(iset < 0 || iset >= snap->num_datasets[itype])
    {
      PyErr_SetString(_read_eagleError, "Dataset index out of range");
      return NULL;  
    }

  /* Return the number of datasets */
  return Py_BuildValue("s", snap->dataset_name[itype]+iset*MAX_NAMELEN);
}


/* Split the selected region between processors */
static PyObject *_read_eagle_split_selection(PyObject *self, PyObject *args)
{
  PyObject      *pysnap;
  EagleSnapshot *snap;
  int NTask, ThisTask;

  /* Get snapshot object and communicator size/rank */
  if (!PyArg_ParseTuple(args, "Oii", &pysnap, &ThisTask, &NTask)) return NULL;
  snap = (EagleSnapshot *) PTR_FROM_OBJECT(pysnap);
  
  /* Call the C function */
  if(split_selection(snap, ThisTask, NTask) != 0)
    {
      PyErr_SetString(_read_eagleError, 
		      "Unable to split selection (already split or bad parameters)");
      return NULL;
    }
  
  /* Success! */
  Py_INCREF(Py_None);
  return Py_None;
}


