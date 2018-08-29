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
find_package(Mosek REQUIRED)
find_package(NLOPT REQUIRED)
find_package(Sophus REQUIRED)


#find_package(NLOPT REQUIRED QUIET)
#include_directories(${NLOPT_INCLUDE_DIRS})

# supereight
set (SE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../../supereight")
message(${SE_MODULE_PATH})

# configure supereight core directoriy
set (SE_CORE_INCLUDE_PATH "${SE_MODULE_PATH}/se_core/include" "${SE_MODULE_PATH}/se_core/include/se")
message(${SE_CORE_INCLUDE_PATH})
# configure supereight denseslam directory
set (SE_DENSESLAM_INCLUDE_PATH "${SE_MODULE_PATH}/se_denseslam/include" "${SE_MODULE_PATH}/se_denseslam/include/se")
message(${SE_DENSESLAM_INCLUDE_PATH})
#set variables
set (SE_INCLUDE_PATH "${SE_CORE_INCLUDE_PATH}" "${SE_DENSESLAM_INCLUDE_PATH}")

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
target_include_directories(${MAV_LIBRARY} PRIVATE ${EIGEN3_INCLUDE_DIR} ${GLOG_INCLUDE_DIRS} ${NLOPT_INCLUDE_DIRS} ${SE_INCLUDE_PATH} ${SOPHUS_INCLUDE_DIR})
target_include_directories(${MAV_LIBRARY} PUBLIC ${MAV_INCLUDE_PATH} ${MOSEK_INCLUDE_DIR} $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${MAV_LIBRARY}>
        $<INSTALL_INTERFACE:include/${MAV_LIBRARY}> glog)
target_link_libraries(${MAV_LIBRARY} ${MOSEK_LIBRARY} ${NLOPT_LIBRARIES})

add_executable(mav_main src/main.cpp)
target_include_directories(mav_main PUBLIC  ${MAV_INCLUDE_PATH} ${EIGEN3_INCLUDE_DIR} ${GLOG_INCLUDE_DIRS} ${NLOPT_INCLUDE_DIRS} ${SE_INCLUDE_PATH} ${SOPHUS_INCLUDE_DIR})
target_link_libraries(mav_main PUBLIC ${MAV_LIBRARY} glog ${NLOPT_LIBRARIES})

#set_target_properties(${MAV_LIBRARY} PROPERTIES PUBLIC_HEADER ${MAV_INCLUDE_PATH})

install(TARGETS ${MAV_LIBRARY} EXPORT ${PROJECT_NAME}Targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

install(DIRECTORY ${MAV_INCLUDE_PATH}/${MAV_LIBRARY} DESTINATION include)