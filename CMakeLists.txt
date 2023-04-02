cmake_minimum_required(VERSION 3.21)

project(quantum_sim LANGUAGES CXX)


add_library(Eigen INTERFACE)
target_include_directories(Eigen INTERFACE externals/eigen-4.3.0)

add_library(ioncpp)

target_compile_features(ioncpp PRIVATE cxx_std_17)
set_target_properties(ioncpp PROPERTIES CXX_STANDARD 17)
set_target_properties(ioncpp PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(ioncpp PUBLIC Eigen)

find_package(OpenMP QUIET)
if(OpenMP_FOUND)
	message("OpenMP Enabled")
	target_link_libraries(Eigen INTERFACE OpenMP::OpenMP_CXX)
endif(OpenMP_FOUND)

add_subdirectory(ioncpp)



add_subdirectory(externals/pybind11)
pybind11_add_module(ionsim ionsim.cpp)
target_compile_features(ionsim PRIVATE cxx_std_17)
set_target_properties(ionsim PROPERTIES CXX_STANDARD 17)
set_target_properties(ionsim PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(ionsim PRIVATE ioncpp)



add_executable(main)

target_sources(main PRIVATE main.cpp)

target_compile_features(main PRIVATE cxx_std_17)
set_target_properties(main PROPERTIES CXX_STANDARD 17)
set_target_properties(main PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(main PRIVATE ioncpp)