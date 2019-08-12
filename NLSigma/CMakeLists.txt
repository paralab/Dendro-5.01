cmake_minimum_required(VERSION 2.8)
project(NLSM)

# options for nlsm
option(NLSM_NONLINEAR "enable nonlinear version" OFF)
option(NLSM_COMPARE_WITH_ANALYTICAL_SOL "compare with the analytical sol. of the linear wave equation" OFF)

if(NLSM_COMPARE_WITH_ANALYTICAL_SOL)
    set(NLSM_NONLINEAR OFF)
    add_definitions(-DNLSM_COMPARE_WITH_ANALYTICAL_SOL)
endif()

if(NLSM_NONLINEAR)
    add_definitions(-DNLSM_NONLINEAR)
endif()


set(NLSM_INC include/parameters.h
        include/rk4nlsm.h
        include/nlsmUtils.h
        include/rhs.h
        include/derivs.h
        include/profile_params.h
        )

set(NLSM_SRC src/rk4nlsm.cpp
        src/parameters.cpp
        src/nlsmUtils.cpp
        src/rhs.cpp
        src/derivs.cpp
        src/profile_params.cpp
        )

set(SOURCE_FILES src/nlsm.cpp  include/nlsm.h ${NLSM_INC} ${NLSM_SRC})
add_executable(nlsmSolver ${SOURCE_FILES})
target_include_directories(nlsmSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(nlsmSolver PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(nlsmSolver PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(nlsmSolver PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(nlsmSolver PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(nlsmSolver PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(nlsmSolver PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(nlsmSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(nlsmSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(nlsmSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(nlsmSolver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(nlsmSolver PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(nlsmSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
