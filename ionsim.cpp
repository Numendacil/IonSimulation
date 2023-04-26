#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include <ioncpp/RungeKutta.hpp>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace ioncpp;
using namespace pybind11::literals;

// Wrapper around the original function, mainly dealing with storage orders
auto func_wrapper(
	const pybind11::EigenDRef<const ArrayType>& init_r,
	const pybind11::EigenDRef<const ArrayType>& init_v,
	const Eigen::Ref<const Eigen::Array<data_t, Eigen::Dynamic, 1>>& charge,
	const Eigen::Ref<const Eigen::Array<data_t, Eigen::Dynamic, 1>>& mass,
	size_t step,
	data_t time_start,
	data_t time_end,
	std::function<
		ArrayType(
			const Eigen::Ref<const ArrayType>& r, 
			const Eigen::Ref<const ArrayType>& v, 
			data_t t
		)
	> force
)
{
	return CalcTrajRK(
		init_r, 
		init_v, 
		charge, 
		mass, 
		step, 
		time_start, 
		time_end, 
		force
	);
}

// NOLINTNEXTLINE
PYBIND11_MODULE(ionsim, m) 
{
	m.def(
		"calculate_trajectory", 
		&func_wrapper,
		"init_r"_a.noconvert(),
		"init_v"_a.noconvert(),
		"charge"_a.noconvert(),
		"mass"_a.noconvert(),
		"step"_a,
		"time_start"_a,
		"time_end"_a,
		"force"_a.noconvert(),
		pybind11::return_value_policy::move
	);

	#ifdef VERSION_INFO
		m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
	#else
		m.attr("__version__") = "dev";
	#endif
}