# DELPHI Nano-AOD
#
# CMakeLists.txt
#
# Create DELPHI NanoAOD
#
# Dietrich Liko <Dietrich.Liko@oeaw.ac.at>

cmake_minimum_required(VERSION 3.26)

project(DELPHI-NanoAOD VERSION 0.1 LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(DELPHI REQUIRED)
find_package(ROOT REQUIRED)

# fetch argparse
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
)
FetchContent_MakeAvailable(yaml-cpp)


add_subdirectory(delphi-analysis)
add_subdirectory(delphi-nanoaod)
add_subdirectory(scripts)

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(delphi-nanoaod stdc++fs)
endif()
