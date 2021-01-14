cmake_minimum_required(VERSION 3.17)
project(DeadCommunity)

set(CMAKE_CXX_STANDARD 14)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug build")
    set(CPLEX_WIN_PLATFORM "x64_windows_msvc14/stat_mdd")
elseif (CMAKE_BUILD_TYPE MATCHES Release)
    message("Release build")
    # OpenMP
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
    set(CPLEX_WIN_PLATFORM "x64_windows_msvc14/stat_mda")
else ()
    message("Some other build type")
endif ()

# CPLEX
set(CPLEX_ROOT_DIR "C:/Program Files/IBM/ILOG/CPLEX_Studio1210")
set(CPLEX_INCLUDE_DIR "${CPLEX_ROOT_DIR}/cplex/include")
set(CPLEX_CONCERT_INCLUDE_DIR "${CPLEX_ROOT_DIR}/concert/include")
set(CPLEX_CP_INCLUDE_DIR "${CPLEX_ROOT_DIR}/cpoptimizer/include")
set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIR} ${CPLEX_CP_INCLUDE_DIR})
include_directories(${CPLEX_INCLUDE_DIRS})

add_definitions("-DIL_STD")

set(CPLEXLIBDIR "${CPLEX_ROOT_DIR}/cplex/lib/${CPLEX_WIN_PLATFORM}")
set(CONCERTLIBDIR "${CPLEX_ROOT_DIR}/concert/lib/${CPLEX_WIN_PLATFORM}")
set(CPLIBDIR "${CPLEX_ROOT_DIR}/cpoptimizer/lib/${CPLEX_WIN_PLATFORM}")
set(CCLNDIRS ${CPLEXLIBDIR} ${CONCERTLIBDIR} ${CPLIBDIR})
link_directories(${CCLNDIRS})

set(CPLEXLIBS ilocplex concert cplex12100 cp)
# END CPLEX

include_directories(.)
include_directories(sfmt)

file(GLOB SOURCES *.cpp sfmt/SFMT.c)

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES})
target_link_libraries(${CMAKE_PROJECT_NAME} ${CPLEXLIBS})
