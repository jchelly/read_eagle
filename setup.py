#
# Script to build the read_eagle python extension module
#
# Run with
#
# python ./setup.py install --prefix=<dir>
#

# Edit this to set the location of HDF5. This directory should
# contain the HDF5 lib/ and include/ directories.
hdf5_location = "/usr/"

from distutils.core import setup, Extension
import numpy.distutils.misc_util

numpy_include_dir = numpy.distutils.misc_util.get_numpy_include_dirs()
idirs             = numpy_include_dir + [hdf5_location+"/include"]
ldirs             = [hdf5_location+"/lib"]
rdirs             = ldirs

read_eagle_module = Extension('_read_eagle',
                              sources = ['./src/_read_eagle.c','./src/read_eagle.c'],
                              libraries=["hdf5"],
                              include_dirs =idirs,
                              library_dirs =ldirs,
                              runtime_library_dirs=rdirs)

setup (name         = 'ReadEagle',
       version      = '1.0',
       description  = 'Code for reading P-H key sorted Eagle snapshots',
       ext_modules  = [read_eagle_module],
       py_modules   = ['read_eagle'],
       package_dir  = {'' : 'src'})
