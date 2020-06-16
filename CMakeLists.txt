# SCRAPP - Species Counting on Reference trees viA Phylogenetic Placements
# Copyright (C) 2020 Pierre Barbera, Lucas Czech and HITS gGmbH
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Pierre Barbera <pierre.barbera@h-its.org>
# Exelixis Lab, Heidelberg Institute for Theoretical Studies
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# adapted from the main CMakeLists of gappa, and the apps-CMakeLists of genesis

# --------------------------------------------------------------------------------------------------
#   CMake Init
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 2.8.7 FATAL_ERROR)

# Fun with colors!
if (NOT WIN32)
    string(ASCII 27 Esc)

    set(ColorBlue   "${Esc}[34m") # Build information.
    set(ColorGreen  "${Esc}[32m") # Optional additions.
    set(ColorYellow "${Esc}[33m") # Tipps for the user.
    set(ColorRed    "${Esc}[31m") # Warnings and Errors.
    set(ColorEnd    "${Esc}[0m")
endif()

set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set (CMAKE_DISABLE_SOURCE_CHANGES  ON)

# --------------------------------------------------------------------------------------------------
#   Project Definition and Requirements
# --------------------------------------------------------------------------------------------------

project ( scrapp-apps CXX )

# --------------------------------------------------------------------------------------------------
#   Build Options
# --------------------------------------------------------------------------------------------------

# The build type is determined by the standard cmake variable CMAKE_BUILD_TYPE.
# Thus, if you want to change the build type, call `cmake .. -DCMAKE_BUILD_TYPE=DEBUG`
IF(CMAKE_BUILD_TYPE STREQUAL "")
    set (CMAKE_BUILD_TYPE RELEASE)
ENDIF()

# Print build type. Useful for debugging user issues.
message (STATUS "${ColorBlue}scrapp-apps build type: ${CMAKE_BUILD_TYPE}${ColorEnd}")

# Store everything in `bin`.
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )

# Option to enable building a fully static version, that is, statical linking of system libraries.
# This works with GCC, but might cause trouble with OpenMP and/or Clang, so it is disabled by
# default. Most users who compile on their own probaby won't need the fully portable version.
# Genesis on the other hand is always linked statically, in order to be able to move around the
# final binary at least on the same system.
option ( APPS_BUILD_STATIC  "Build the static versions of the executables." OFF )

if( APPS_BUILD_STATIC AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    string(ASCII 27 Esc)
    message( STATUS "${Esc}[31m"
        "Compiler is not GCC. This might cause trouble with statically linking system libraries "
        "on Linux systems. If the compilation fails with some linking error, try a different compiler, "
        "or try to not build statically by deactivating the CMake option 'APPS_BUILD_STATIC'. "
    "${Esc}[0m" )
endif()

# --------------------------------------------------------------------------------------------------
#   Compiler and Linker Options
# --------------------------------------------------------------------------------------------------

# Set the warn flags to a very high level - except unknown pragmas, as this causes needless
# warnings with OpenMP and other pragma-based techniques.
set (WARN_FLAGS              "-Wall -Wextra -Wno-unknown-pragmas -pedantic -pedantic-errors")
set (CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++11 ${WARN_FLAGS}")

set (CMAKE_CXX_FLAGS_DEBUG   "-O2 -DDEBUG -g -ggdb3 -D_GLIBCXX_DEBUG" )
set (CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG" )

IF(APPS_BUILD_STATIC)
    # Use different ways of telling CMake to build static, just to be sure.
    SET(BUILD_SHARED_LIBRARIES OFF)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    link_libraries("-static")
ENDIF()
message( STATUS "Static linking of system libraries: ${APPS_BUILD_STATIC}")

# --------------------------------------------------------------------------------------------------
#   Dependencies Settings
# --------------------------------------------------------------------------------------------------

# Add Genesis as dependency.
add_subdirectory(${PROJECT_SOURCE_DIR}/programs/genesis)

# Use everything that Genesis exports, just to be sure that we use the same setup.
add_definitions( ${GENESIS_DEFINITIONS} )
include_directories( ${GENESIS_INCLUDE_DIR} )
set( CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          ${GENESIS_C_FLAGS}")
set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${GENESIS_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GENESIS_EXE_LINKER_FLAGS}" )

# Genesis exports the OpenMP variables. We use it to check and give advice in case it is not found.
if( NOT OPENMP_FOUND )
    string(ASCII 27 Esc)
    message( STATUS "${Esc}[31m"
        "OpenMP was not found or has been deactivated. This results in a considerably slower "
        "runtime for the program. Try to use a different compiler. Alternatively, if you are sure "
        "your compiler supports OpenMP, you can add the needed compiler flags to the CMake script."
        # " See BUILD.md for details."
    "${Esc}[0m" )
endif()

# Set rpath for all targets, so that they work out of the box (libgenesis.so in parent directory),
# as well as when relocating with the lib being in the same directory.
# Also, we need to set two versions of ORIGIN here (one and two dollar signs).
# This is because on some systems it seems to be escaped, on some not...
# See https://cmake.org/pipermail/cmake/2008-January/019290.html
# and https://cmake.org/Wiki/CMake_RPATH_handling for details.
# Tipp: Use `ldd binary` and `readelf -d binary` to check RPATH settings.
set( CMAKE_INSTALL_RPATH "$ORIGIN/..:$$ORIGIN/..:$ORIGIN:$$ORIGIN" )
set( CMAKE_BUILD_WITH_INSTALL_RPATH TRUE )

# Add dependencies that have CMake build
# include_directories( "libs/CLI11/include" )


# ------------------------------------------------------------------------------
#   Sources
# ------------------------------------------------------------------------------

# Find all source files in this directory and compile them into individual binaries, each
# linked to the shared lib. All files should thus contain their own main function.
file( GLOB apps_sources ${PROJECT_SOURCE_DIR}/src/*.cpp )

include_directories( ${PROJECT_SOURCE_DIR}/deps/genesis/libs )

# ------------------------------------------------------------------------------
#   Build Applications
# ------------------------------------------------------------------------------

if( apps_sources )

    message( STATUS "${ColorBlue}Building applications:${ColorEnd}" )
    foreach( app_src ${apps_sources} )
        get_filename_component( app_name ${app_src} NAME_WE )
        message( STATUS "  ${ColorBlue}${app_name}${ColorEnd}" )

        add_executable( ${app_name} ${app_src} )

        target_link_libraries( ${app_name} ${GENESIS_LINK_LIBRARIES} )

        # Link against any external libraries, e.g. Pthreads.
        target_link_libraries (${app_name} ${GENESIS_INTERNAL_LINK_LIBRARIES})

    endforeach()

else()

    message( STATUS "No scapp applications applications found" )

endif()
