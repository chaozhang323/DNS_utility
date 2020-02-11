cmake_minimum_required(VERSION 3.5.1)
project(test_kdtree)
enable_language(Fortran)


## library
## HDF5
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS Fortran)
include_directories(${HDF5_INCLUDE_DIRS})

## LAPACK
find_package(LAPACK)

## include
#include_directories(include/)

## compiler
set(CMAKE_Fortran_COMPILER gfortran)

## compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message("Enabled GNU Fortran.")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffree-line-length-none")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -cpp")
    message("Enabled Intel Fortran.")
endif()

#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fcheck=all -Wall")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DVERBOSE")

message("${HDF5_Fortran_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES}")
message("Compiling with flags ${CMAKE_Fortran_FLAGS}")
## compile
file(GLOB SOURCES ../*.f90 *.f90 )
add_executable(test_kdtree ${SOURCES})
target_link_libraries(test_kdtree PUBLIC ${HDF5_Fortran_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES})
