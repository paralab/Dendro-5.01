cmake_minimum_required(VERSION 2.8)
project(QUADGRAV)

# options for quadgrav

set(QUADGRAV_INC include/parameters.h
        include/rk4quadgrav.h
        include/quadgravUtils.h
        include/rhs.h
        include/derivs.h
        include/profile_params.h
        )

set(QUADGRAV_SRC src/rk4quadgrav.cpp
        src/parameters.cpp
        src/quadgravUtils.cpp
        src/rhs.cpp
        src/derivs.cpp
        src/profile_params.cpp
        )

set(SOURCE_FILES src/quadgrav.cpp  include/quadgrav.h ${QUADGRAV_INC} ${QUADGRAV_SRC})
add_executable(quadgravSolver ${SOURCE_FILES})
install(TARGETS quadgravSolver RUNTIME DESTINATION bin)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(quadgravSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(quadgravSolver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(quadgravSolver PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(quadgravSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
