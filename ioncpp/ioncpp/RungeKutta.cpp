#include "RungeKutta.hpp"
#include <iostream>

namespace ioncpp
{

using namespace std;

pair<vector<ArrayType>, vector<ArrayType>> CalcTrajRK(
	const Eigen::Ref<const ArrayType>& init_r,
	const Eigen::Ref<const ArrayType>& init_v,
	size_t step,
	data_t time_start,
	data_t time_end,
	AcclCallback accl
)
{
	assert(time_end >= time_start);	// NOLINT(*-pro-bounds-array-to-pointer-decay)
	data_t dt = (time_end - time_start) / (data_t)step;

	ArrayType a(DIM, init_r.cols());
	vector<ArrayType> r_ret, v_ret;
	r_ret.reserve(step);
	v_ret.reserve(step);

	ArrayType r = init_r;
	ArrayType v = init_v;

	ArrayType r_tmp(DIM, init_r.cols());
	ArrayType v_tmp(DIM, init_r.cols());

	ArrayType r_k1(DIM, init_r.cols());
	ArrayType r_k2(DIM, init_r.cols());
	ArrayType r_k3(DIM, init_r.cols());
	ArrayType r_k4(DIM, init_r.cols());

	ArrayType v_k1(DIM, init_r.cols());
	ArrayType v_k2(DIM, init_r.cols());
	ArrayType v_k3(DIM, init_r.cols());
	ArrayType v_k4(DIM, init_r.cols());

	for (size_t i = 0; i < step; i++)
	{
		double t = time_start + dt * (double)step;
		
		v_k1 = accl(r, v, t) * dt;
		r_k1 = v * dt;
		r_tmp = r + r_k1 / 2.0;
		v_tmp = v + v_k1 / 2.0;

		v_k2 = accl(r_tmp, v_tmp, t + dt / 2.0) * dt;
		r_k2 = v_tmp * dt;
		r_tmp = r + r_k2 / 2.0;
		v_tmp = v + v_k2 / 2.0;

		v_k3 = accl(r_tmp, v_tmp, t + dt / 2.0) * dt;
		r_k3 = v_tmp * dt;
		r_tmp = r + r_k3 / 2.0;
		v_tmp = v + v_k3 / 2.0;

		v_k4 = accl(r_tmp, v_tmp, t + dt / 2.0) * dt;
		r_k4 = v_tmp * dt;
		v += v_k1 / 6.0 + v_k2 / 3.0 + v_k3 / 3.0 + v_k4 / 6.0;
		r += r_k1 / 6.0 + r_k2 / 3.0 + r_k3 / 3.0 + r_k4 / 6.0;

		r_ret.push_back(r);
		v_ret.push_back(v);
	}
	return {r_ret, v_ret};
}


}