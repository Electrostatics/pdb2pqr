cmake_minimum_required (VERSION 3.13)

#set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_FIND_FRAMEWORK NEVER)

find_package(Python3 COMPONENTS Interpreter Development)
if(${Python3_FOUND})
    message(STATUS "***** Python3 FOUND")
    message(STATUS "***** Python3 include path is: ${Python3_INCLUDE_DIRS}")
    message(STATUS "***** Python3 library path is: ${Python3_LIBRARIES}")
    message(STATUS "***** Python3 library dir  is: ${Python3_LIBRARY_DIRS}")
endif()

find_program(VIRTUALENV virtualenv)
if(NOT VIRTUALENV)
    message(FATAL_ERROR "Could not find `virtualenv` in PATH")
endif()
set(VIRTUALENV ${VIRTUALENV} -p ${Python3_EXECUTABLE})

# Generate the virtualenv and ensure it's up to date.
add_custom_command(
    OUTPUT venv
    COMMAND ${VIRTUALENV} venv
)
add_custom_command(
    OUTPUT venv.stamp
    DEPENDS venv requirements.txt
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt requirements.txt
    COMMAND ./venv/bin/pip install -r requirements.txt --upgrade
)

#
#  Now that we have setup the virtual environment with NumPy, we want to find the
#  the NumPy include files using the new virtual environment
#
set(Python3_ROOT_DIR ${CMAKE_BINARY_DIR}/pdb2pqr/venv)
find_package(Python3 COMPONENTS Interpreter Development NumPy)
if(${Python3_FOUND})
    message(STATUS "***** VENV Python3 FOUND")
    message(STATUS "***** VENV Python3 include path is: ${Python3_INCLUDE_DIRS}")
    message(STATUS "***** VENV Python3 library path is: ${Python3_LIBRARIES}")
    message(STATUS "***** VENV Python3 library dir  is: ${Python3_LIBRARY_DIRS}")
    message(STATUS "***** VENV Python3 NumPy   path is: ${Python3_NumPy_INCLUDE_DIRS}")
endif()

add_library(Algorithms
    ${CMAKE_CURRENT_SOURCE_DIR}/pdb2pqr/pdb2pka/substruct/Algorithms.cpp
)
target_include_directories(Algorithms PUBLIC ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(Algorithms ${Python3_LIBRARIES})
#Python3_add_library(Algorithms MODULE ${CMAKE_CURRENT_SOURCE_DIR}/pdb2pqr/pdb2pka/substruct/Algorithms.cpp)

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
