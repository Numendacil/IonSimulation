#ifndef ION_SIMULATION_RUNGE_KUTTA_HPP
#define ION_SIMULATION_RUNGE_KUTTA_HPP

#include <Eigen/Dense>
#include <vector>

namespace ioncpp
{

constexpr std::size_t DIM = 3;

using data_t = double;

using ArrayType = Eigen::Array<data_t, Eigen::Dynamic, DIM>;
using ArrayTypeRef = Eigen::Ref<const ArrayType>;

using ForceCallback = std::function<ArrayType(const ArrayTypeRef& r, const ArrayTypeRef& v, data_t t)>;

/**
 * @brief Calculate ion trajectory using 4th Runge-Kutta
 * 
 * @param init_r the initial position of ions
 * @param init_v the initial velocity of ions
 * @param charge the charge of ions
 * @param mass the mass of ions
 * @param step iteration step
 * @param time_start start time
 * @param time_end end time
 * @param force dynamic function for external force based on input position, velocity, time, excluding
 * the Coulomb interaction between ions which has already taken into account
 * @return std::pair<std::vector<ArrayType>, std::vector<ArrayType>>
 */
std::pair<std::vector<ArrayType>, std::vector<ArrayType>> CalcTrajRK(
	const ArrayTypeRef& init_r,
	const ArrayTypeRef& init_v,
	const Eigen::Ref<const Eigen::Array<data_t, Eigen::Dynamic, 1>>& charge,
	const Eigen::Ref<const Eigen::Array<data_t, Eigen::Dynamic, 1>>& mass,
	size_t step,
	data_t time_start,
	data_t time_end,
	ForceCallback force
);

}

#endif