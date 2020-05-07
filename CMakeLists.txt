cmake_minimum_required(VERSION 3.12)

project(pdb2pqr)

include(CMakePrintHelpers)

# set(CMAKE_VERBOSE_MAKEFILE ON)

###############################################################################
###  CMAKE_FIND_FRAMEWORK by default looks in MACOSX Frameworks FIRST whenever
###  find_package is used. This allows cmake to find other version of packages
###  if they are installed (e.g. /usr/local or $HOME/.anaconda3)
###############################################################################
set(CMAKE_FIND_FRAMEWORK NEVER)

###############################################################################
###  We need to find python3 so that we can install a virtual environment
###  with all the packages required.
###############################################################################
#find_package(Python3 COMPONENTS Interpreter REQUIRED)
#if(${Python3_FOUND})
#    message(STATUS "***** Python3 FOUND          : ${Python3_EXECUTABLE}")
#endif()

###############################################################################
###  Set up a python3 virtual environment so we don't pollute the user's
###  default python environment and so their environment doesn't interfere
###  with pdb2pqr
###############################################################################
#find_program(VIRTUALENV virtualenv)
#if(NOT VIRTUALENV)
#    message(FATAL_ERROR "Could not find `virtualenv` in PATH")
#endif()

###############################################################################
###  Generate the virtualenv and ensure it's up to date.
###############################################################################
#if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/venv")
#    execute_process(
#        COMMAND ${VIRTUALENV} -p ${Python3_EXECUTABLE} venv
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        RESULT_VARIABLE VIRTUALENV_RESULT)
#    execute_process(
#        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt requirements.txt
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        RESULT_VARIABLE COPY_RESULT)
#    set(ENV{VIRTURAL_ENV} ${CMAKE_CURRENT_BINARY_DIR}/venv)
#    execute_process(
#        COMMAND pip3 install -r requirements.txt --upgrade
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        RESULT_VARIABLE PIP_RESULT)
#endif()

#
#  Now that we have setup the virtual environment with NumPy, we want to find the
#  the NumPy include files using the new virtual environment
#
unset(Python3_FOUND) 
set(Python3_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/venv")
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
if(${Python3_FOUND})
    message(STATUS "***** VENV Python3 FOUND          : ${Python3_EXECUTABLE}")
    message(STATUS "***** VENV Python3 include path is: ${Python3_INCLUDE_DIRS}")
    message(STATUS "***** VENV Python3 library path is: ${Python3_LIBRARIES}")
    message(STATUS "***** VENV Python3 library dir  is: ${Python3_LIBRARY_DIRS}")
    message(STATUS "***** VENV Python3 NumPy   path is: ${Python3_NumPy_INCLUDE_DIRS}")
endif()

# TODO: This is only for debugging
return()
message(STATUS "***** THIS SHOULD NOT HAPPEN")

# Build command line to run py.test.
set(PYTEST
    ${CMAKE_CURRENT_BINARY_DIR}/venv/bin/python
    ${CMAKE_CURRENT_BINARY_DIR}/venv/bin/py.test
)

add_custom_target(Tests ALL
    DEPENDS venv.stamp
    SOURCES requirements.txt
)

###############################################################################
#####  Create Tests
###############################################################################
if (${ENABLE_TESTS})
    message(STATUS "** Tests will be built")
    enable_testing()

    foreach(testsourcefile
        ${CMAKE_CURRENT_SOURCE_DIR}/pytest/00_pytest_test.py
    )
        get_filename_component(t_filename ${testsourcefile} NAME_WE)
        add_test(NAME ${t_filename}_tests
            COMMAND ${PYTEST} ${testsourcefile}
        )
    endforeach(testsourcefile)
endif()

include(CTest)
