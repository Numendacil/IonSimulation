#include <fstream>
#include <iostream>
#include <ion_sim/Runge_Kutta.hpp>

int main()
{
	using namespace ion_sim;
	ArrayType r0(3, 1);
	ArrayType v0(3, 1);
	r0 << 5, 0, 0;
	v0 << 0, 5, 0;

	AcclCallback f = [](const auto& r, const auto& /*v*/, data_t /*t*/) -> ArrayType
	{
		return -r * (25.0 / 5.0);
	};

	auto [r_list, v_list] = runge_kutta(r0, v0, 2000, 0, 20, std::move(f));
	std::ofstream r_file("r_list.txt");
	std::ofstream v_file("v_list.txt");

	Eigen::IOFormat format(Eigen::StreamPrecision, Eigen::DontAlignCols, ",");
	for (const auto& p : r_list)
		r_file << p.transpose().format(format) << std::endl;
	for (const auto& p : v_list)
		v_file << p.transpose().format(format) << std::endl;

	return 0;
}