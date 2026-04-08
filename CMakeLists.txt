cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# project name and enable C++ support
project(nlodis CXX C)

set(CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Git module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_DESCRIBE)
git_local_changes(GIT_LOCAL_CHANGES)

find_package(GSL REQUIRED)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(DEBUG)
endif()

# Cuba library distributed as a part of this project
include_directories("src/Cuba-4.2.2/")

include(CTest)

add_subdirectory(src)

# Several tests read gbw.dat using a relative path from the current working directory.
# Keep a copy in the build root so running tests from build/ works reliably.
configure_file("${CMAKE_SOURCE_DIR}/src/tests/gbw.dat" "${CMAKE_BINARY_DIR}/gbw.dat" COPYONLY)

# Copy the BK datafile to the build directory as well, since the example code uses that
file(COPY "${CMAKE_SOURCE_DIR}/nlobkdatafiles"
    DESTINATION "${CMAKE_BINARY_DIR}")

set_target_properties(nlodis_tests PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )

if (BUILD_TESTING)
    add_test(NAME unit_tests COMMAND $<TARGET_FILE:nlodis_tests>)
    set_tests_properties(unit_tests PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
endif()

set_target_properties(nlodis PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )

set_target_properties(simple_example PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin" )