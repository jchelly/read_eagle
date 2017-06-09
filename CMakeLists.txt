cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(read_eagle)
enable_language(C)
enable_language(Fortran)

#
# Override cmake's default behaviour of stripping library
# paths from executables
#
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#
# Find HDF5
#
find_package(HDF5 COMPONENTS C REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_LIBRARY_DIRS})

#
# Determine how to call C from Fortran
#
include(FortranCInterface)
FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/include/FC.h MACRO_NAMESPACE "FC_")

#
# Where to output executables, libraries etc
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include/)

#
# Subdirectories with the code
#
add_subdirectory(src)
add_subdirectory(examples)

#
# Packaging stuff
# This allows you to generate a source tar file with "make package_source"
#
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_PACKAGE_FILE_NAME "read_eagle" CACHE INTERNAL "tarball basename")
set(CPACK_SOURCE_IGNORE_FILES
"~$"
"\\\\.gitignore$"
"\\\\/.git/"
)
include(CPack)

