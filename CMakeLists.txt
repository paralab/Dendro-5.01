cmake_minimum_required(VERSION 2.8)
project(Dendro-5.01)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# For now we just make it compulsory to have LAPACK installed.
#Later we will make it possible if LAPACK is not present to automaticall install before compiling dendro5





if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()




option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off" ON)
option(HILBERT_ORDERING "Use Hilbert ordering instead of Morton" ON)
option(RUN_WEAK_SCALING "Run Entire Weak Scaling" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" ON)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" ON)
option(DIM_2 "use the two dimentional sorting" OFF)
option(PROFILE_TREE_SORT "Profile the tree sort code, construction and balancing. " OFF)
option(WITH_BLAS_LAPACK "build using BLAS and LAPACk" OFF)
option(MANUAL_BLAS_LAPACK "configure BLAS and LAPACK Manually" OFF)
#option(KWAY "K parameter for alltoallv_kway" 128)
set(KWAY 128 CACHE INT 128)
set(NUM_NPES_THRESHOLD 2 CACHE INT 2)
#message(${KWAY})


if(WITH_BLAS_LAPACK)
    add_definitions(-DWITH_BLAS_LAPACK)
if(MANUAL_BLAS_LAPACK)
    if("$ENV{BLAS}" STREQUAL "")
        message("Environment Variable BLAS is not set. Please set it to BLAS directory")
    endif()

    if( "$ENV{LAPACK}" STREQUAL "" )
        message("Enviroment Variable LAPACK is note set. Please set it to LAPACK directory. ")
    endif()
    set(LAPACKE_DIR $ENV{LAPACK}/LAPACKE)
    set(BLAS_LIBS $ENV{BLAS}/lib )
    #set(CBLAS_LIBS $ENV{CBLAS}/lib )
    set(LAPACK_LIBS $ENV{LAPACK}/lib)
    include_directories(${LAPACKE_DIR}/include)
    #set(LAPACKE_DIR /home/milinda/Programs/lapack-3.7.0/LAPACKE)
    #set(BLAS_LIBS /home/milinda/Programs/BLAS-3.7.0/lib )
    #set(CBLAS_LIBS $ENV{CBLAS}/lib )
    #set(LAPACK_LIBS /home/milinda/Programs/lapack-3.7.0/lib)

    set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
    set(LAPACK_LIBRARIES ${LAPACK_LIBS}/liblapacke.a ${LAPACK_LIBS}/liblapack.a ${BLAS_LIBS}/libblas.a libgfortran.a libquadmath.a)
    set(LINK_FLAGS "${LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
else ()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(LAPACKE_DIR $ENV{LAPACK}/LAPACKE)
    include_directories(${LAPACK_INCLUDE_DIR})
    #link_directories(${MKL_LIBRARIES})
    set(LINK_FLAGS "${LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
endif()

endif()

if(DIM_2)
    add_definitions(-DDIM_2)
endif()

if(PROFILE_TREE_SORT)
    add_definitions(-DPROFILE_TREE_SORT)
endif()


if(USE_64BIT_INDICES)
    add_definitions(-DUSE_64BIT_INDICES)
    #message('Configuring 64BIT indices')
endif()

if(HILBERT_ORDERING)
    add_definitions(-DHILBERT_ORDERING)
    #message (Configuring with HIlbert)
endif()

if(RUN_WEAK_SCALING)
    add_definitions(-DRUN_WEAK_SCALING)
    #message('Configuring with Weak Scaling')
endif()

if(ALLTOALLV_FIX)
    add_definitions(-DALLTOALLV_FIX)
    add_definitions(-DKWAY=${KWAY})
endif()

if(SPLITTER_SELECTION_FIX)
    add_definitions(-DSPLITTER_SELECTION_FIX)
    add_definitions(-DNUM_NPES_THRESHOLD=${NUM_NPES_THRESHOLD})
endif()

set(LAPACK_INC /home/milinda/Programs/lapack-3.7.0/LAPACKE/include/)
include_directories(${PROJECT_BINARY_DIR}
        include
        FEM/include
        ODE/include
        IO/vtk/include
        include/test
        ${MPI_INCLUDE_PATH}
        examples/include
        )



add_library(hcurve include/hcurvedata.h src/hcurvedata.cpp)
add_library(genPts include/genPts_par.h src/genPts_par.cpp)
add_library(treeNode include/TreeNode.h src/TreeNode.cpp)
add_library(treeNode2vtk include/treenode2vtk.h src/treenode2vtk.cpp)

set(DENDRO_INCLUDE_FILES
        include/point.h
        include/seqUtils.h
        include/seqUtils.tcc
        include/ompUtils.h
        include/ompUtils.tcc
        include/dendro.h
        include/binUtils.h
        include/hcurvedata.h
        include/genPts_par.h
        include/TreeNode.h
        include/treenode2vtk.h
        include/parUtils.h
        include/parUtils.tcc
        include/dendro.h
        include/colors.h
        include/dtypes.h
        include/sfcSort.h
        include/sfcSearch.h
        include/test/testUtils.h
        include/test/testUtils.tcc
        include/octUtils.h
        include/dendroIO.h
        include/mathUtils.h
        include/mathUtils.tcc
        include/stencil.h
        include/dollar.hpp
        )
set(DENDRO_SOURCE_FILES
        src/point.cpp
        src/binUtils.cpp
        src/hcurvedata.cpp
        src/genPts_par.cpp
        src/TreeNode.cpp
        src/treenode2vtk.cpp
        src/parUtils.cpp
        src/octUtils.cpp
        src/dendroIO.cpp
        )


set(IO_INC IO/vtk/include/oct2vtk.h
           IO/vtk/include/json.hpp
        )

set(IO_SRC IO/vtk/src/oct2vtk.cpp
        )


add_library(dendro5 ${DENDRO_INCLUDE_FILES} ${DENDRO_SOURCE_FILES} ${IO_INC} ${IO_SRC})

set(SOURCE_FILES Hilbert/src/genHTables.cpp include/radix.h)
add_executable(genHtables ${SOURCE_FILES})


set(SOURCE_FILES examples/src/ptsSort.cpp )
add_executable(ptsSort ${SOURCE_FILES})
target_link_libraries(ptsSort dendro5 ${MPI_LIBRARIES} m)

set(SOURCE_FILES
        examples/include/testOctUtils.h
        examples/src/testOCtUtils.cpp
        )
add_executable(testOctUtils ${SOURCE_FILES})
target_link_libraries(testOctUtils dendro5 ${MPI_LIBRARIES} m)


set(SOURCE_FILES
        examples/include/octUtilsBenchmark.h
        examples/src/octUtilsBenchmark.cpp
        )


add_executable(octUtilBenchmark ${SOURCE_FILES})
target_link_libraries(octUtilBenchmark dendro5 ${MPI_LIBRARIES} m)


set(SOURCE_FILES
        examples/include/bucketBench.h
        examples/src/bucketBench.cpp
        )


add_executable(bucketBench ${SOURCE_FILES})
target_link_libraries(bucketBench dendro5 ${MPI_LIBRARIES} m)

