#
# Script to build the read_eagle python extension module
#
# Run with
#
# python ./setup.py install --prefix=<dir>
#
# Edit the line below to set the location of HDF5. This directory
# should contain the HDF5 lib/ and include/ directories. If the HDF5
# bin directory is in your $PATH then you can use
#
# h5cc -showconfig
#
# to determine where HDF5 is installed.
#
hdf5_location = "/usr/"

import sys
from distutils.core import setup, Extension
import numpy.distutils.misc_util

numpy_include_dir = numpy.distutils.misc_util.get_numpy_include_dirs()
idirs             = numpy_include_dir + [hdf5_location+"/include"]
ldirs             = [hdf5_location+"/lib"]

if sys.platform.startswith("win"):
    # On Windows, must not specify run time search path
    rdirs = []
    # For static HDF5 library
    #extra_compile_args = ["-D_HDF5USEDLL_",]
    # For dynamic HDF5 library
    extra_compile_args = ["-DH5_BUILT_AS_DYNAMIC_LIB",]
else:
    # Set runtime library search path on non-Windows systems
    rdirs = ldirs
    # No need for extra args in this case
    extra_compile_args = []

read_eagle_module = Extension('_read_eagle',
                              sources = ['./src/_read_eagle.c','./src/read_eagle.c'],
                              libraries=["hdf5"],
                              include_dirs =idirs,
                              library_dirs =ldirs,
                              runtime_library_dirs=rdirs,
                              extra_compile_args=extra_compile_args,
                          )

setup (name         = 'ReadEagle',
       version      = '1.0',
       description  = 'Code for reading P-H key sorted Eagle snapshots',
       ext_modules  = [read_eagle_module],
       py_modules   = ['read_eagle'],
       package_dir  = {'' : 'src'})
