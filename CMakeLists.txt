cmake_minimum_required(VERSION 2.8)

add_definitions("/MP")

project(shadow-fraction)
SET(BOOST_ROOT C:/installer/boost_1_67_0/)
find_package(Boost REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/include" "${Boost_INCLUDE_DIRS}")
link_directories("${PROJECT_SOURCE_DIR}/lib" "${Boost_LIBRARY_DIRS}")

add_subdirectory(src)
