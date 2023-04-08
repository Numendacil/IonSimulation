#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include <ioncpp/RungeKutta.hpp>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

// NOLINTNEXTLINE
PYBIND11_MODULE(ionsim, m) 
{
	m.def("calculate_trajectory", &ioncpp::CalcTrajRK);

	#ifdef VERSION_INFO
		m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
	#else
		m.attr("__version__") = "dev";
	#endif
}