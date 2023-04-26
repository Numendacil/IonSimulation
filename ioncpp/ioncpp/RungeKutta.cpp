#include "RungeKutta.hpp"

#include <iostream>

namespace ioncpp
{

using namespace std;

namespace
{

ArrayType ColoumbInteraction(
	const ArrayTypeRef& r, 
	const Eigen::Ref<const Eigen::Array<data_t, Eigen::Dynamic, 1>>& charge
)
{
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
			/ dist2.array().pow(1.5) * (charge.matrix() * charge.transpose().matrix()).array()
		).rowwise().sum();
	}
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
		
		v_k1 = (force(r, v, t)+ ColoumbInteraction(r, charge)).colwise() / mass * dt;
		r_k1 = v * dt;
		r_tmp = r + r_k1 / 2.0;
		v_tmp = v + v_k1 / 2.0;

		v_k2 = (force(r_tmp, v_tmp, t + dt / 2.0)+ ColoumbInteraction(r, charge)).colwise() / mass * dt;
		r_k2 = v_tmp * dt;
		r_tmp = r + r_k2 / 2.0;
		v_tmp = v + v_k2 / 2.0;

		v_k3 = (force(r_tmp, v_tmp, t + dt / 2.0)+ ColoumbInteraction(r, charge)).colwise() / mass * dt;
		r_k3 = v_tmp * dt;
		r_tmp = r + r_k3 / 2.0;
		v_tmp = v + v_k3 / 2.0;

		v_k4 = (force(r_tmp, v_tmp, t + dt)+ ColoumbInteraction(r, charge)).colwise() / mass * dt;
		r_k4 = v_tmp * dt;
		v += v_k1 / 6.0 + v_k2 / 3.0 + v_k3 / 3.0 + v_k4 / 6.0;
		r += r_k1 / 6.0 + r_k2 / 3.0 + r_k3 / 3.0 + r_k4 / 6.0;

		r_ret.push_back(r);
		v_ret.push_back(v);
	}
	return {r_ret, v_ret};
}


}