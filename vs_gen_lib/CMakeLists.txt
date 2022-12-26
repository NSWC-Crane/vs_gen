# ----------------------------------------------------------------------------
#  Root CMake file for the virtual scene generation library
# ----------------------------------------------------------------------------
#
# To perform an out of source build do the following on the command line:
# 1. navigate to the root project directory (in this case DfD2)
# 2. mkdir build
# 3. cd build
# 4. depending on your system (windows/linux) type the following
# Example CMake build line for windows & VS2015
# windows: cmake -G "Visual Studio 14 2015 Win64" -DUSE_AVX_INSTRUCTIONS=ON ..
# linux: cmake -DUSE_AVX_INSTRUCTIONS=ON ..
#
# 5. To compile:
# windows: cmake --build . --config Release
# linux: cmake --build . --config Release -- -j10
#

cmake_minimum_required(VERSION 3.9)

project(vs_gen_lib VERSION 1.0 DESCRIPTION "Shared library for generating virtual scenes for DfD training using opencv")
set(CMAKE_CXX_STANDARD 11)

set(LIB_NAME "vs_gen")
set(default_build_type "Release")

set(PLATFORM "$ENV{PLATFORM}" CACHE INTERNAL "Copied from environment variable")
if("${PLATFORM}" STREQUAL "")
    set(PLATFORM "NOT_SET")
    message("The PLATFORM Environment variable is not set on this system.  Setting it to:")
endif()
message("Platform: " ${PLATFORM})

# find the root of the main project
string(FIND ${PROJECT_SOURCE_DIR} "/" INDEX REVERSE)
string(SUBSTRING ${PROJECT_SOURCE_DIR} 0 ${INDEX} PR_ROOT)

# Build Options
option(USE_AVX_INSTRUCTIONS  "Compile your program with AVX instructions"  ON)

# ----------------------------------------------------------------------------
if(WIN32)
    if(${PLATFORM} STREQUAL "Laptop")
        set(COMMON_ROOT "D:/Common")
        set(TURB_ROOT "D:/Projects/turbulence_sim")
        set(RYML_ROOT "D:/Projects/rapidyaml")
    elseif(${PLATFORM} STREQUAL "LaptopN")
        set(COMMON_ROOT "C:/Projects/Common")
        set(TURB_ROOT "C:/Projects/turbulence_sim")
        set(RYML_ROOT "C:/Projects/rapidyaml")
    else()
        set(COMMON_ROOT "C:/Projects/Common")
        set(TURB_ROOT "C:/Projects/turbulence_sim")
        set(RYML_ROOT "C:/Projects/rapidyaml")
    endif()

    set(RYML_LIB "${RYML_ROOT}/build/Release/ryml.lib")

else()
    set(HOME "$ENV{HOME}" CACHE INTERNAL "Copied from environment variable")

    if(${PLATFORM} STREQUAL "HPC")
        set(COMMON_ROOT "${HOME}/Projects/Common")
        set(TURB_ROOT "${HOME}/Projects/turbulence_sim")
        set(RYML_ROOT "${HOME}/Projects/rapidyaml")
    elseif(${PLATFORM} STREQUAL "MainGear")
        set(COMMON_ROOT "${HOME}/Projects/Common")
        set(TURB_ROOT "${HOME}/Projects/turbulence_sim")
        set(RYML_ROOT "${HOME}/Projects/rapidyaml")
    else()
        set(COMMON_ROOT "${HOME}/Projects/Common")
        set(TURB_ROOT "${HOME}/Projects/turbulence_sim")
        set(RYML_ROOT "${HOME}/Projects/rapidyaml")
    endif()

    set(RYML_LIB "${RYML_ROOT}/build/libryml.so")

endif()

include(${COMMON_ROOT}/cmake_utils/use_sse_avx.cmake)


# ----------------------------------------------------------------------------
# get the source files for the Library
file(GLOB SOURCES "src/*.cpp")

set(project_includes "include"
    "${PR_ROOT}/common/include"
    "${TURB_ROOT}/common/include"
    "${COMMON_ROOT}/include"
    "${RYML_ROOT}/include"
    )

set(project_libs ${project_libs} ${RYML_LIB})

# ----------------------------------------------------------------------------
# find OpenCV
find_package(OpenCV REQUIRED)
if (OpenCV_FOUND)
    set(project_includes ${project_includes} ${OpenCV_INCLUDE_DIRS})
    set(project_libs ${project_libs} ${OpenCV_LIBS})
else()
    message("OpenCV not found, so we won't build the project.")
endif()


# ----------------------------------------------------------------------------
message(" ")
message("--------------------------------------------------------------------------------")
message("Source Files: ")
foreach(S ${SOURCES})
    message("  " ${S})
endforeach()
message("--------------------------------------------------------------------------------")

message(" ")

message("--------------------------------------------------------------------------------")
message("Project Includes: ")
foreach(S ${project_includes})
    message("  " ${S})
endforeach()
message("--------------------------------------------------------------------------------")

message(" ")

message("--------------------------------------------------------------------------------")
message("Project Libraries: ")
foreach(S ${project_libs})
    message("  " ${S})
endforeach()
message("--------------------------------------------------------------------------------")

message(" ")

# ----------------------------------------------------------------------------
message("Building Library...")

include_directories(${project_includes})

if(WIN32)

    add_definitions(-DLIB_EXPORTS)

    add_library(${LIB_NAME} SHARED ${SOURCES})
    set_target_properties(${LIB_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
    target_link_libraries(${LIB_NAME} ${project_libs})

    if(MSVC AND (${CMAKE_VERSION} VERSION_GREATER "3.6.0"))
        message("MSVC Version: " ${MSVC_VERSION})
        target_compile_options(${LIB_NAME} PUBLIC "/MP" "/bigobj")
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${LIB_NAME})
    endif() 

elseif(UNIX)

    add_definitions(-DLIB_EXPORTS)
    # Change these to point to the right location for a particular machine
    #add_definitions("-fPIC")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    #find_package(X11 REQUIRED)
    add_library(${LIB_NAME} SHARED ${SOURCES})
    set_target_properties(${LIB_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
    #target_link_libraries(${LIB_NAME} ${X11_LIBRARIES})
    target_link_libraries(${LIB_NAME} ${project_libs})

elseif(APPLE)
    message(FATAL_ERROR "Apple -> Not Even!")

endif()

# ----------------------------------------------------------------------------
message(" ")
message(STATUS "Build Complete!")
message(" ")

