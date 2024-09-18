cmake_minimum_required(VERSION 3.12)

# set the project name
project(diagarpack LANGUAGES Fortran)
include(FindBLAS)
include(FindLAPACK)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_Fortran_FLAGS "-DBLOCK_CUTOFF -Ofast -g -Wall")

find_package(BLAS)
find_package(LAPACK)

find_library(ARPACK_LIBRARIES NAMES arpack PATHS /usr/local/lib /opt/homebrew/lib)
if(ARPACK_LIBRARIES)
	MESSAGE(STATUS "Found ARPACK at ${ARPACK_LIBRARIES}")
else(ARPACK_LIBRARIES)
	MESSAGE(FATAL_ERROR "ARPACK not found")
endif(ARPACK_LIBRARIES)

# add the executable
add_executable(rtb2 rtb.f90 diag_arpack.F rtb_mod.F)
target_link_libraries(rtb2 ${ARPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
