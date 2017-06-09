# Put C headers in /include relative to build dir. 
file(COPY read_eagle.h DESTINATION ${CMAKE_BINARY_DIR}/include/)
include_directories(${CMAKE_BINARY_DIR}/include/)

# Create the shared library
add_library (read_eagle SHARED read_eagle.c read_eagle_f.c read_eagle_fortran.f90)

# Link to HDF5
target_link_libraries(read_eagle ${HDF5_LIBRARIES})

# Make sure that the library and Fortran module get installed if we do "make install"
install (TARGETS read_eagle DESTINATION lib)
install (DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
