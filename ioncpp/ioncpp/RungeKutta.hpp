#ifndef ION_SIMULATION_RUNGE_KUTTA_HPP
#define ION_SIMULATION_RUNGE_KUTTA_HPP

#include <Eigen/Dense>
#include <vector>

namespace ioncpp
{

constexpr std::size_t DIM = 3;

using data_t = double;

using ArrayType = Eigen::Array<data_t, DIM, Eigen::Dynamic, Eigen::RowMajor>;
using AcclCallback = std::function<ArrayType(const Eigen::Ref<const ArrayType>& r, const Eigen::Ref<const ArrayType>& v, data_t t)>;

/**
 * @brief Calculate ion trajectory using 4th Runge-Kutta
 * 
 * @param init_r the initial position of ions
 * @param init_v the initial velocity of ions
 * @param step iteration step
 * @param time_start start time
 * @param time_end end time
 * @param accl dynamic function for acceleration, input position, velocity, time
 * @return std::pair<std::vector<ArrayType>, std::vector<ArrayType>>
 */
std::pair<std::vector<ArrayType>, std::vector<ArrayType>> CalcTrajRK(
	const Eigen::Ref<const ArrayType>& init_r,
	const Eigen::Ref<const ArrayType>& init_v,
	size_t step,
	data_t time_start,
	data_t time_end,
	AcclCallback accl
);

}

#endif