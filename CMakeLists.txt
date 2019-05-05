cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PROJECT_FAZIA)

#----------------------------------------------------------------------------

set(useROFAZIA true)
set(FAZIA_SOURCE_DIRECTORIES src)
set(FAZIA_EXECUTABLES example)
set(FAZIA_SCRIPTS example1.mac)

#----------------------------------------------------------------------------

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

if(useROFAZIA)
  list(APPEND FAZIA_SOURCE_DIRECTORIES)
  list(APPEND FAZIA_EXECUTABLES)
endif(useROFAZIA)

set(FAZIA_SOURCES_DIRECTORIES_FULL)
foreach(_directory ${FAZIA_SOURCE_DIRECTORIES})
  list(APPEND FAZIA_SOURCE_DIRECTORIES_FULL ${PROJECT_SOURCE_DIR}/${_directory})
endforeach(_directory)

set(FAZIA_SOURCES)
foreach(_directory ${FAZIA_SOURCE_DIRECTORIES_FULL})
  list(APPEND FAZIA_SOURCES ${_directory}/*.cc)
endforeach(_directory)

set(FAZIA_HEADERS)
foreach(_directory ${FAZIA_SOURCE_DIRECTORIES_FULL})
  list(APPEND FAZIA_HEADERS ${_directory}/*.hh)
endforeach(_directory)

include(${Geant4_USE_FILE})
include_directories(${FAZIA_SOURCE_DIRECTORIES_FULL})

if(useROFAZIA)
  EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROFAZIA_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROFAZIA_CXX_FLAGS}")
endif(useROFAZIA)

file(GLOB sources ${FAZIA_SOURCES})
file(GLOB headers ${FAZIA_HEADERS})

foreach(_executable ${FAZIA_EXECUTABLES})
  add_executable(${_executable} ${_executable}.cc ${sources} ${headers})
if(useROFAZIA)
  EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROFAZIA_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_EXE_LINKER_FLAGS ${ROFAZIA_LD_FLAGS})
endif(useROFAZIA)
  target_link_libraries(${_executable} ${Geant4_LIBRARIES})
endforeach(_executable)

foreach(_script ${FAZIA_SCRIPTS})
  configure_file( ${PROJECT_SOURCE_DIR}/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach(_script)

add_custom_target(PROJECT_FAZIA DEPENDS ${FAZIA_EXECUTABLES})

foreach(_executable ${FAZIA_EXECUTABLES})
  install(TARGETS ${_executable} DESTINATION bin)
endforeach(_executable)
