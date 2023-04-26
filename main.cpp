#include <fstream>
#include <iostream>
#include <ioncpp/RungeKutta.hpp>

int main()
{
	using namespace ioncpp;
	Eigen::Array<data_t, Eigen::Dynamic, DIM> r0(2, 3);
	Eigen::Array<data_t, Eigen::Dynamic, DIM> v0(2, 3);
	Eigen::Array<data_t, Eigen::Dynamic, 1> q(2, 1);
	Eigen::Array<data_t, Eigen::Dynamic, 1> m(2, 1);
	r0 << 1, 0, 0,
		-1, 0, 0;
	v0 << 0, 1, 0,
		0, -1, 0;
	q << 1, -1;
	m << 1, 1;

	ForceCallback f = [&m](const auto& r, const auto& /*v*/, data_t /*t*/) -> ArrayType
	{
		return -((r * 0.75).colwise() * m);
	};

	auto [r_list, v_list] = CalcTrajRK(
		r0, v0, 
		q,
		m,
		20000, 
		0, 20, 
		std::move(f));
	std::ofstream r_file("r_list.txt");
	std::ofstream v_file("v_list.txt");

	Eigen::IOFormat format(Eigen::StreamPrecision, Eigen::DontAlignCols, ",");
	for (const auto& p : r_list)
		r_file << p.transpose().format(format) << std::endl;
	for (const auto& p : v_list)
		v_file << p.transpose().format(format) << std::endl;

	return 0;
}