cmake_minimum_required(VERSION 3.16)
project(rbac_maintenance_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

# EvalMaxSAT2022 with IPAMIR API (incremental MaxSAT solver with assumptions)
set(EVALMAXSAT2022_DIR "${CMAKE_SOURCE_DIR}/third_party/EvalMaxSAT2022")

if(EXISTS "${EVALMAXSAT2022_DIR}")
    message(STATUS "EvalMaxSAT2022 found at ${EVALMAXSAT2022_DIR}")

    set(EVALMAXSAT2022_INCLUDES
        ${CMAKE_SOURCE_DIR}/third_party
        ${EVALMAXSAT2022_DIR}
        ${EVALMAXSAT2022_DIR}/lib/EvalMaxSAT/src
        ${EVALMAXSAT2022_DIR}/lib/MaLib/src
        ${EVALMAXSAT2022_DIR}/lib/cadical/src
    )

    set(EVALMAXSAT2022_LIB "${EVALMAXSAT2022_DIR}/libipamirEvalMaxSAT2022.a")

    add_library(rbac_maintenance
        src/SATSolver.cpp
        src/BMF.cpp
        src/BMFLocalSearch.cpp
    )

    target_include_directories(rbac_maintenance PUBLIC
        src
        ${EVALMAXSAT2022_INCLUDES}
    )

    if(EXISTS "${EVALMAXSAT2022_LIB}")
        target_link_libraries(rbac_maintenance ${EVALMAXSAT2022_LIB})
        message(STATUS "EvalMaxSAT2022 library: ${EVALMAXSAT2022_LIB}")
    else()
        message(FATAL_ERROR "libipamirEvalMaxSAT2022.a not found. Run ./build.sh first to compile EvalMaxSAT2022")
    endif()

    find_package(ZLIB REQUIRED)
    find_package(Threads REQUIRED)
    target_link_libraries(rbac_maintenance ZLIB::ZLIB Threads::Threads)

else()
    message(FATAL_ERROR "EvalMaxSAT2022 not found at ${EVALMAXSAT2022_DIR}")
endif()

# Main executable
add_executable(rbac_main src/main.cpp)
target_link_libraries(rbac_main rbac_maintenance)
