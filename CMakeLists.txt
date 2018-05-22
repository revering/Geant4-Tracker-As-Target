cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4DM)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
   find_package(Geant4 REQUIRED ui_all vis_all)
else()
   find_package(Geant4 REQUIRED)
endif()

find_package(GSL REQUIRED)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJEXT_SOURCE_DIR}/include/*.hh)

add_executable(G4DM G4DM.cc ${sources} ${headers})
target_link_libraries(G4DM ${Geant4_LIBRARIES} ${GSL_LIBRARIES})

set(G4DM_SCRIPTS
#  G4DM.in
#  G4DM.out
  init_vis.mac
#  run1.mac
#  run2.mac
  vis.mac
  )

foreach(_script ${G4DM_SCRIPTS})
   configure_file(
      ${PROJECT_SOURCE_DIR}/${_script}
      ${PROJECT_BINARY_DIR}/${_script}
      COPYONLY
      )
endforeach()


install(TARGETS G4DM DESTINATION bin)
