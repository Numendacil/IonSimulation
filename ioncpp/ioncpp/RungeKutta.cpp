#include "RungeKutta.hpp"

#include <iostream>
#include <chrono>

namespace ioncpp
{

using namespace std;
using namespace std::literals;

namespace
{

// chrono::microseconds elapsed1 = 0us;

ArrayType CoulombInteraction(
	const ArrayTypeRef& r, 
	const Eigen::Ref<const Eigen::Array<data_t, Eigen::Dynamic, 1>>& charge
)
{
	// auto begin = chrono::steady_clock::now();

	const long N = r.rows();
	Eigen::Matrix<data_t, Eigen::Dynamic, Eigen::Dynamic> dist2 = 
	(
		r.rowwise().squaredNorm().matrix() 
		* Eigen::Matrix<data_t, 1, Eigen::Dynamic>::Ones(1, N)
	)
	+ 
	(
		Eigen::Matrix<data_t, Eigen::Dynamic, 1>::Ones(N, 1) 
		* r.rowwise().squaredNorm().transpose().matrix()
	)
	- 2 * (r.matrix() * r.matrix().transpose());
	dist2.diagonal() = Eigen::Vector<data_t, Eigen::Dynamic>::Ones(N);

	ArrayType result(N, 3);
	for (uint8_t i = 0; i < 3; i++)
	{
		result.col(i) = (
			(r.col(i).matrix() * Eigen::Matrix<data_t, 1, Eigen::Dynamic>::Ones(1, N)
			- Eigen::Matrix<data_t, Eigen::Dynamic, 1>::Ones(N, 1) * r.col(i).transpose().matrix()).array()
			/ (dist2.array().sqrt() * dist2.array()) * (charge.matrix() * charge.transpose().matrix()).array()
		).rowwise().sum();
	}

	// auto end = std::chrono::steady_clock::now();
	// elapsed1 += std::chrono::duration_cast<std::chrono::microseconds>(end - begin);

	return result;
}

}

pair<vector<ArrayType>, vector<ArrayType>> CalcTrajRK(
	const ArrayTypeRef& init_r,
	const ArrayTypeRef& init_v,
	const Eigen::Ref<const Eigen::Array<data_t, Eigen::Dynamic, 1>>& charge,
	const Eigen::Ref<const Eigen::Array<data_t, Eigen::Dynamic, 1>>& mass,
	size_t step,
	data_t time_start,
	data_t time_end,
	ForceCallback force
)
{
	// auto begin = chrono::steady_clock::now();
	// auto tmp = chrono::steady_clock::now();
	// chrono::microseconds elapsed2 = 0us;

	assert(time_end >= time_start);	// NOLINT(*-pro-bounds-array-to-pointer-decay)
	data_t dt = (time_end - time_start) / (data_t)step;

	ArrayType a(init_r.rows(), DIM);
	vector<ArrayType> r_ret, v_ret;
	r_ret.reserve(step);
	v_ret.reserve(step);

	ArrayType r = init_r;
	ArrayType v = init_v;

	ArrayType r_tmp(init_r.rows(), DIM);
	ArrayType v_tmp(init_r.rows(), DIM);

	ArrayType r_k1(init_r.rows(), DIM);
	ArrayType r_k2(init_r.rows(), DIM);
	ArrayType r_k3(init_r.rows(), DIM);
	ArrayType r_k4(init_r.rows(), DIM);

	ArrayType v_k1(init_r.rows(), DIM);
	ArrayType v_k2(init_r.rows(), DIM);
	ArrayType v_k3(init_r.rows(), DIM);
	ArrayType v_k4(init_r.rows(), DIM);

	for (size_t i = 0; i < step; i++)
	{
		double t = time_start + dt * (double)i;
		
		// tmp = std::chrono::steady_clock::now();
		v_k1 = (force(r, v, t)+ CoulombInteraction(r, charge)).colwise() / mass * dt;
		// elapsed2 += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - tmp);
		r_k1 = v * dt;
		r_tmp = r + r_k1 / 2.0;
		v_tmp = v + v_k1 / 2.0;

		// tmp = std::chrono::steady_clock::now();
		v_k2 = (force(r_tmp, v_tmp, t + dt / 2.0)+ CoulombInteraction(r_tmp, charge)).colwise() / mass * dt;
		// elapsed2 += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - tmp);
		r_k2 = v_tmp * dt;
		r_tmp = r + r_k2 / 2.0;
		v_tmp = v + v_k2 / 2.0;

		// tmp = std::chrono::steady_clock::now();
		v_k3 = (force(r_tmp, v_tmp, t + dt / 2.0)+ CoulombInteraction(r_tmp, charge)).colwise() / mass * dt;
		// elapsed2 += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - tmp);
		r_k3 = v_tmp * dt;
		r_tmp = r + r_k3;
		v_tmp = v + v_k3;

		// tmp = std::chrono::steady_clock::now();
		v_k4 = (force(r_tmp, v_tmp, t + dt)+ CoulombInteraction(r_tmp, charge)).colwise() / mass * dt;
		// elapsed2 += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - tmp);
		r_k4 = v_tmp * dt;
		v += v_k1 / 6.0 + v_k2 / 3.0 + v_k3 / 3.0 + v_k4 / 6.0;
		r += r_k1 / 6.0 + r_k2 / 3.0 + r_k3 / 3.0 + r_k4 / 6.0;

		r_ret.push_back(r);
		v_ret.push_back(v);
	}

	// auto end = chrono::steady_clock::now();

	// std::cout << "Time elapsed in CoulombInteraction: " << elapsed1.count() << "[µs]" << std::endl;
	// std::cout << "Time elapsed in force callback: " << elapsed2.count() - elapsed1.count() << "[µs]" << std::endl;
	// std::cout << "Time elapsed in CalcTrajRK: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

	return {r_ret, v_ret};
}


}