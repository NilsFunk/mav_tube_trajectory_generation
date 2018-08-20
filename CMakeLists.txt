cmake_minimum_required(VERSION 3.1)

project(mav_tube_trajectory_generation)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

add_compile_options(-std=c++14 -march=native)

set (MAV_LIBRARY "${PROJECT_NAME}")
set (MAV_MODULE_PATH "${PROJECT_SOURCE_DIR}")
set (MAV_SRC_PATH  "${MAV_MODULE_PATH}/src" )
set (MAV_INCLUDE_PATH  "${MAV_MODULE_PATH}/include")


# Build se_denseslam lib
find_package(Eigen3 REQUIRED)
find_package(Glog REQUIRED QUIET)


#find_package(NLOPT REQUIRED QUIET)
#include_directories(${NLOPT_INCLUDE_DIRS})

#############
# LIBRARIES #
#############
add_library(${MAV_LIBRARY} SHARED
  src/motion_defines.cpp
  src/polynomial.cpp
  src/segment.cpp
  src/timing.cpp
  src/trajectory.cpp
  src/vertex.cpp
  src/rpoly/rpoly_ak1.cpp
)

#target_link_libraries(${PROJECT_NAME} )
target_include_directories(${MAV_LIBRARY} PRIVATE ${EIGEN3_INCLUDE_DIR} ${GLOG_INCLUDE_DIRS})
target_include_directories(${MAV_LIBRARY} PUBLIC ${MAV_INCLUDE_PATH}   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${MAV_LIBRARY}>
        $<INSTALL_INTERFACE:include/${MAV_LIBRARY}>)

#set_target_properties(${MAV_LIBRARY} PROPERTIES PUBLIC_HEADER ${MAV_INCLUDE_PATH})

install(TARGETS ${MAV_LIBRARY} EXPORT ${PROJECT_NAME}Targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

install(DIRECTORY ${MAV_INCLUDE_PATH}/${MAV_LIBRARY} DESTINATION include)