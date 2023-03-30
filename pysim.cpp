#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include <ion_sim/Runge_Kutta.hpp>

// NOLINTNEXTLINE
PYBIND11_MODULE(pysim, m) 
{
	m.def("runge_kutta", &ion_sim::runge_kutta);
}