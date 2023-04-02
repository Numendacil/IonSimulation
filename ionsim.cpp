#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include <ioncpp/RungeKutta.hpp>

// NOLINTNEXTLINE
PYBIND11_MODULE(ionsim, m) 
{
	m.def("calculate_trajectory", &ioncpp::CalcTrajRK);

	m.attr("__version__") = "dev";
}