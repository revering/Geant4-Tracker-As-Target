cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TTarg)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
   find_package(Geant4 REQUIRED ui_all vis_all)
else()
   find_package(Geant4 REQUIRED)
endif()

find_package(GSL REQUIRED )

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJEXT_SOURCE_DIR}/include/*.hh)

add_executable(TTarg tracker_as_target.cc ${sources} ${headers})
target_link_libraries(TTarg ${Geant4_LIBRARIES} ${GSL_LIBRARIES})

set(TTarg_SCRIPTS
  TTarg.in
  TTarg.out
  init_vis.mac
#  run1.mac
#  run2.mac
  vis.mac
  )

foreach(_script ${TTarg_SCRIPTS})
   configure_file(
      ${PROJECT_SOURCE_DIR}/${_script}
      ${PROJECT_BINARY_DIR}/${_script}
      COPYONLY
      )
endforeach()


install(TARGETS TTarg DESTINATION bin)
