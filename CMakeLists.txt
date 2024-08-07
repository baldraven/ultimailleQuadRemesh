#cmake_minimum_required(VERSION 3.12...3.26)
#
#project( playground LANGUAGES C CXX )
#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#
#get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
#if(NOT is_multi_config AND NOT (CMAKE_BUILD_TYPE OR DEFINED ENV{CMAKE_BUILD_TYPE}))
#  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release default")
#endif()
#
#find_package(OpenMP COMPONENTS CXX)
#
#if (MSVC)
#    # warning level 4 (and all warnings as errors, /WX)
#    add_compile_options(/W4)
#else()
#    # lots of warnings and all warnings as errors
#    add_compile_options(-Wall -Wextra -pedantic)
#endif()

## Clone Ultimaille
#include(FetchContent)
#FetchContent_Declare(
#  ultimaille
#  GIT_REPOSITORY https://github.com/ssloy/ultimaille
#  GIT_TAG        v1.1.0
#)
#FetchContent_MakeAvailable(ultimaille)
#include_directories(${ultimaille_SOURCE_DIR})

#add_subdirectory(src)



cmake_minimum_required(VERSION 3.12...3.26)

project(ultimailleBunin LANGUAGES C CXX )
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)

message("OSName: " ${OSName})
message("System name: " ${CMAKE_SYSTEM_NAME})


# Get config
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT is_multi_config AND NOT (CMAKE_BUILD_TYPE OR DEFINED ENV{CMAKE_BUILD_TYPE}))
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release default")
endif()

# Get development version of ultimaille directly from git repository 
include(FetchContent)
FetchContent_Declare(
  ultimaille
  GIT_REPOSITORY https://github.com/ssloy/ultimaille
  # GIT_TAG        v1.1.0
  GIT_TAG        geom
)
FetchContent_MakeAvailable(ultimaille)

FetchContent_Declare(
  param_parser
  GIT_REPOSITORY https://github.com/ultimaille/param-parser
  GIT_TAG        v1.0.0
)
FetchContent_MakeAvailable(param_parser)



# Add ultimaille directory to get ultimaille CMakeLists
include_directories(${ultimaille_SOURCE_DIR})
include_directories(${param_parser_SOURCE_DIR})

# Add executable binary
add_executable(main src/main.cpp)

# Link binary with ultimaille / param-parser
target_link_libraries(main  ${CMAKE_DL_LIBS} param-parser ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)

IF (NOT WIN32)
	target_link_libraries(main m)
ENDIF()

if (NOT DEFINED OSName)
  set(OSName ${CMAKE_SYSTEM_NAME})
endif()

set_target_properties(
  main
    PROPERTIES 
        OUTPUT_NAME "main_${OSName}"
)
