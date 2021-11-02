# preset that will explicitly request gcc/g++ compilers with support for MPI and OpenMP

set(CMAKE_CXX_COMPILER "g++" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "gcc" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O2 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(MPI_CXX "g++" CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
set(MPI_C "gcc" CACHE STRING "" FORCE)
set(MPI_C_COMPILER "mpicc" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-Wall -g" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O2 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(MPI_Fortran "gfortran" CACHE STRING "" FORCE)
set(MPI_Fortran_COMPILER "mpifort" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG "-Wall -g -std=f2003" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O2 -DNDEBUG -std=f2003" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -std=f2003" CACHE STRING "" FORCE)
unset(HAVE_OMP_H_INCLUDE CACHE)

set(OpenMP_C "gcc" CACHE STRING "" FORCE)
set(OpenMP_C_FLAGS "-fopenmp" CACHE STRING "" FORCE)
set(OpenMP_C_LIB_NAMES "gomp" CACHE STRING "" FORCE)
set(OpenMP_CXX "g++" CACHE STRING "" FORCE)
set(OpenMP_CXX_FLAGS "-fopenmp" CACHE STRING "" FORCE)
set(OpenMP_CXX_LIB_NAMES "gomp" CACHE STRING "" FORCE)
set(OpenMP_omp_LIBRARY "libgomp.so" CACHE PATH "" FORCE)
