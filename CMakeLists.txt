# ##############################################################################
# CMAKE CONFIGURATION FILE FOR DENDRO DIRECTORY
#
# This file in particular contains all of the options that apply specificially
# to the Dendro directory as it is compiled. This includes the example solvers
# as well as the library itself.
#
# NOTE: there is a CMakeLists-NonLibraryDefs.txt file which contains more options.
# The options in there must be included in the root directory of your project
# that you are using this library with. This will configure everything for building,
# but the options must be propagated to your own solvers and projects.
#
# ##############################################################################

cmake_minimum_required(VERSION 3.15)
project(Dendro-5.0-Library C CXX Fortran)


# ==========================
# DENDRO LIBRARY OPTIONS
# ==========================
option(USE_64BIT_INDICES "Use 64-Bit indices (Reverts to 32-bit if OFF)" ON)
option(HILBERT_ORDERING "Use Hilbert ordering instead of Morton" ON)
option(RUN_WEAK_SCALING "Run Entire Weak Scaling" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" OFF)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" OFF)
option(DIM_2 "Use the two dimensional sorting" OFF)
option(PROFILE_TREE_SORT "Profile the tree sort code" OFF)
option(WITH_BLAS_LAPACK "Build using BLAS and LAPACK" ON)
option(DENDRO_VTK_BINARY "Write VTK/VTU files in binary mode" ON)
option(DENDRO_VTK_ZLIB_COMPRES "Compress VTK/VTU (requires Binary)" OFF)
option(ALLTOALL_SPARSE "Use isend/irecv for ghost exchange" OFF)
option(ENABLE_DENDRO_PROFILE_COUNTERS "Enable dendro internal profile counters" OFF)
option(RK_SOLVER_OVERLAP_COMM_AND_COMP "Enable RK solver comm/comp overlap" ON)
option(WITH_CUDA "Build dendro with CUDA" OFF)
option(BUILD_WITH_PETSC "Build dendro with PETSc" OFF)
option(USE_FD_INTERP_FOR_UNZIP "Use FD style interpolation for unzip" OFF)
option(MPIX_CUDA_AWARE_SUPPORT "Allow MPI to be aware of CUDA support" OFF)
option(DVEC_ZERO_ALLOC "Initialize dvec memory to zero" OFF)

# Configuration Variables
set(KWAY 128 CACHE STRING "K parameter for alltoallv_kway")
set(OCT2BLK_COARSEST_LEV 0 CACHE STRING "Coarsest Level")
set(NUM_NPES_THRESHOLD 2 CACHE STRING "NPES Threshold")



# ==========================
# DEPENDENCIES
# ==========================
include(FetchContent)

# OpenMP
find_package(OpenMP REQUIRED)

# MPI
find_package(MPI REQUIRED)

# GSL
find_package(GSL REQUIRED)

# spdlog allows for logging!
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.14.1
)
FetchContent_MakeAvailable(spdlog)

# make sure we have BhaHAHA included
if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/BHaHAHA")
  add_subdirectory(lib/BHaHAHA)
endif()


# ==========================
# SOURCE FILES
# ==========================
message(STATUS "Configuring Dendro Library...")

set(DENDRO_SRC
    src/dendro.cpp src/point.cpp src/binUtils.cpp src/hcurvedata.cpp
    src/genPts_par.cpp src/TreeNode.cpp src/treenode2vtk.cpp src/parUtils.cpp
    src/mesh.cpp src/octUtils.cpp src/dendroIO.cpp src/block.cpp
    src/sc_flops.cpp src/profiler.cpp src/dendroProfileParams.cpp src/oda.cpp
    src/odaUtils.cpp src/sub_oda.cpp src/meshUtils.cpp src/waveletRefEl.cpp
    test/src/meshTestUtils.cpp src/aeh_bhahaha.cpp src/logger.cpp
)

set(FEM_SRC
    FEM/examples/src/laplace.cpp FEM/src/basis.cpp FEM/src/refel.cpp
    FEM/src/tensor.cpp FEM/src/workspace.cpp)

set(ODE_SRC
    ODE/src/rk.cpp ODE/src/subSM.cpp ODE/src/rkMaxwell.cpp
    ODE/src/rkTransport.cpp ODE/src/enutsOp.cpp)

set(IO_SRC IO/vtk/src/oct2vtk.cpp IO/vtk/src/checkPoint.cpp)

set(ZLIB_SRC
    IO/zlib/src/adler32.c IO/zlib/src/compress.c IO/zlib/src/crc32.c
    IO/zlib/src/deflate.c IO/zlib/src/gzclose.c IO/zlib/src/gzlib.c
    IO/zlib/src/gzread.c IO/zlib/src/gzwrite.c IO/zlib/src/infback.c
    IO/zlib/src/inffast.c IO/zlib/src/inflate.c IO/zlib/src/inftrees.c
    IO/zlib/src/trees.c IO/zlib/src/uncompr.c IO/zlib/src/zutil.c)

if(WITH_CUDA)
  list(APPEND DENDRO_SRC src/cudaUtils.cpp src/oda.cpp FEM/examples/src/heatMat.cpp)
endif()

FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.14.1
    )
FetchContent_MakeAvailable(spdlog)

# ==========================
# DENDRO LIBRARY DEFINITION
# ==========================
add_library(dendro5 ${DENDRO_SRC} ${FEM_SRC} ${ODE_SRC} ${IO_SRC} ${ZLIB_SRC})

# add definitions to COMPILE for dendro5
macro(dendro_add_def OPTION_NAME DEF_NAME)
  if(${OPTION_NAME})
    target_compile_definitions(dendro5 PUBLIC ${DEF_NAME})
  endif()
endmacro()

target_compile_definitions(dendro5 PUBLIC -DMATVEC_PROFILE)
target_compile_definitions(dendro5 PUBLIC -DOCT2BLK_COARSEST_LEV=${OCT2BLK_COARSEST_LEV})

dendro_add_def(USE_64BIT_INDICES USE_64BIT_INDICES)
dendro_add_def(HILBERT_ORDERING HILBERT_ORDERING)
dendro_add_def(RUN_WEAK_SCALING RUN_WEAK_SCALING)
dendro_add_def(DIM_2 DIM_2)
dendro_add_def(PROFILE_TREE_SORT PROFILE_TREE_SORT)
dendro_add_def(ALLTOALL_SPARSE ALLTOALL_SPARSE)
dendro_add_def(RK_SOLVER_OVERLAP_COMM_AND_COMP RK_SOLVER_OVERLAP_COMM_AND_COMP)
dendro_add_def(USE_FD_INTERP_FOR_UNZIP USE_FD_INTERP_FOR_UNZIP)
dendro_add_def(DVEC_ZERO_ALLOC DVEC_ZERO_ALLOC)

if(ALLTOALLV_FIX)
  target_compile_definitions(dendro5 PUBLIC ALLTOALLV_FIX KWAY=${KWAY})
endif()

if(SPLITTER_SELECTION_FIX)
  target_compile_definitions(dendro5 PUBLIC SPLITTER_SELECTION_FIX NUM_NPES_THRESHOLD=${NUM_NPES_THRESHOLD})
endif()

if(ENABLE_DENDRO_PROFILE_COUNTERS)
  target_compile_definitions(dendro5 PUBLIC ENABLE_DENDRO_PROFILE_COUNTERS __PROFILE_CTX__ __PROFILE_ETS__)
endif()

# VTK Logic
if(DENDRO_VTK_BINARY)
  target_compile_definitions(dendro5 PUBLIC DENDRO_VTU_BINARY)
  if(DENDRO_VTK_ZLIB_COMPRES)
    target_compile_definitions(dendro5 PUBLIC DENDRO_VTU_ZLIB)
  endif()
else()
  target_compile_definitions(dendro5 PUBLIC DENDRO_VTU_ASCII)
endif()

target_include_directories(dendro5 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/test/include
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/include
    ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include
    ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include
    ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/include
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include
    ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc
)

target_link_libraries(dendro5 PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX GSL::gsl spdlog::spdlog BHaHAHA m)

if(WITH_BLAS_LAPACK)
  target_compile_definitions(dendro5 PUBLIC WITH_BLAS_LAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(dendro5 PUBLIC BLAS::BLAS LAPACK::LAPACK)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    target_compile_options(dendro5 PUBLIC "-mkl")
    target_link_options(dendro5 PUBLIC "-mkl")
  endif()
endif()

if(BUILD_WITH_PETSC)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")
  find_package(PETSc REQUIRED)
  target_compile_definitions(dendro5 PUBLIC BUILD_WITH_PETSC)
  target_include_directories(dendro5 PUBLIC ${PETSC_INCLUDES})
  target_link_libraries(dendro5 PUBLIC ${PETSC_LIBRARIES})
endif()

if(WITH_CUDA)
  set(OCT2BLK_COARSEST_LEV 31 CACHE STRING "Coarsest Level for CUDA" FORCE)
  target_compile_definitions(dendro5 PUBLIC WITH_CUDA OCT2BLK_COARSEST_LEV=31)

  find_package(CUDA REQUIRED)
  target_include_directories(dendro5 PUBLIC ${CUDA_INCLUDE_DIRS})

  if(MPIX_CUDA_AWARE_SUPPORT)
    target_compile_definitions(dendro5 PUBLIC MPIX_CUDA_AWARE_SUPPORT)
  endif()

  # Create the separate CUDA target
  set(CUDA_SRCS
        ${DENDRO_SRC} ${FEM_SRC} ${ODE_SRC} ${IO_SRC} ${ZLIB_SRC}
        GPU/include/derivs_cu.cuh GPU/src/refel_const.cu GPU/src/mesh_gpu.cu)

  cuda_add_library(dendro5_cuda ${CUDA_SRCS})

  # We must manually re-apply properties to the cuda library
  target_include_directories(dendro5_cuda PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/GPU/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include
        ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include
        ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include
        ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc
        ${CUDA_INCLUDE_DIRS}
    )
  # Define flags
  target_compile_definitions(dendro5_cuda PUBLIC WITH_CUDA)
  target_link_libraries(dendro5_cuda PUBLIC ${CUDA_LIBRARIES} MPI::MPI_CXX)
endif()

# ==========================
# Examples and Tests!
# ==========================

option(BUILD_DENDRO_EXAMPLES "Build Dendro Examples" ON)

if(BUILD_DENDRO_EXAMPLES)
  function(dendro_add_example TARGET_NAME)
    add_executable(${TARGET_NAME} ${ARGN})
    # NOTE: linking dendro5 includes **all** flags and includes!
    target_link_libraries(${TARGET_NAME} PRIVATE dendro5 MPI::MPI_CXX m)
  endfunction()

  dendro_add_example(genHtables Hilbert/src/genHTables.cpp include/radix.h)
  dendro_add_example(genSFCMatVecTables Hilbert/src/sfcMatvecTable.cpp)

  dendro_add_example(ptsSort examples/src/ptsSort.cpp)
  dendro_add_example(testOctUtils examples/src/testOCtUtils.cpp)
  dendro_add_example(tstTreeSearch examples/src/tstTreeSearch.cpp)
  dendro_add_example(octUtilBenchmark examples/src/octUtilsBenchmark.cpp)
  dendro_add_example(bucketBench examples/src/bucketBench.cpp)
  dendro_add_example(mfreeBench examples/src/mfreeTest.cpp)
  dendro_add_example(meshBench examples/src/meshBenchmark.cpp)
  dendro_add_example(meshE2NCheck examples/src/meshE2NCheck.cpp)
  dendro_add_example(run_all_tests test/src/run_all_tests.cpp)
  dendro_add_example(benchUtils test/src/benchUtils.cpp)

  dendro_add_example(testLTSBlkVec ODE/examples/testLTSBlkVec.cpp)
  dendro_add_example(dgVecTest FEM/examples/src/dgVecTest.cpp)
  dendro_add_example(testInterp FEM/examples/src/testInterp.cpp)
  dendro_add_example(meshLoop examples/src/meshLooping.cpp)
  dendro_add_example(refEl FEM/examples/src/refelExample.cpp)
  dendro_add_example(refElOctree FEM/examples/src/refelOctreeExample.cpp)
  dendro_add_example(transportEq ODE/examples/transportEq.cpp)
  dendro_add_example(laplaceEq FEM/examples/src/laplace.cpp)
  dendro_add_example(heatEq FEM/examples/src/heatMat.cpp FEM/examples/src/heatVec.cpp FEM/examples/src/heatEq.cpp)
  dendro_add_example(subDAExample FEM/examples/src/subdaExample.cpp)
  dendro_add_example(sfcLaplaceEq FEM/examples/src/sfcLaplace.cpp)

  if(WITH_CUDA)
    # cuda GPU tests
    cuda_add_executable(run_meshgpu_tests ${CMAKE_CURRENT_SOURCE_DIR}/GPU/src/run_meshgpu_tests.cu)
    target_link_libraries(run_meshgpu_tests dendro5 dendro5_cuda)
    target_include_directories(run_meshgpu_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/GPU/include)
  endif()
endif()
