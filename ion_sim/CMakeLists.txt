target_sources(
	ion_sim PRIVATE

	ion_sim/Runge_Kutta.hpp
	ion_sim/Runge_Kutta.cpp
)

target_include_directories(ion_sim PUBLIC .)