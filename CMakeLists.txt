cmake_minimum_required(VERSION 3.21)

project(quantum_sim LANGUAGES CXX)


add_library(Eigen INTERFACE)
target_include_directories(Eigen INTERFACE externals/eigen-4.3.0)

add_library(ion_sim)

target_compile_features(ion_sim PRIVATE cxx_std_17)
set_target_properties(ion_sim PROPERTIES CXX_STANDARD 17)
set_target_properties(ion_sim PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(ion_sim PUBLIC Eigen)

find_package(OpenMP QUIET)
if(OpenMP_FOUND)
	message("OpenMP Enabled")
	target_link_libraries(Eigen INTERFACE OpenMP::OpenMP_CXX)
endif(OpenMP_FOUND)

add_subdirectory(ion_sim)



add_executable(main)

target_sources(main PRIVATE main.cpp)

target_compile_features(main PRIVATE cxx_std_17)
set_target_properties(main PROPERTIES CXX_STANDARD 17)
set_target_properties(main PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(main PRIVATE ion_sim)

add_subdirectory(externals/pybind11)
pybind11_add_module(pysim pysim.cpp)
target_compile_features(pysim PRIVATE cxx_std_17)
set_target_properties(pysim PROPERTIES CXX_STANDARD 17)
set_target_properties(pysim PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(pysim PRIVATE ion_sim)