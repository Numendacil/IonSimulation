import ionsim 
import numpy as np

r0 = np.array([
	[1, 0, 0],
	[-1, 0, 0],
], dtype=np.float64)
v0 = np.array([
	[0, 1, 0],
	[0, -1, 0],
], dtype=np.float64)

charge = np.array([1, -1], dtype=np.float64)
mass = np.array([1, 1], dtype=np.float64)

def force(r: np.ndarray, v: np.ndarray, t: float) -> np.ndarray:
	return - r * 0.75 * mass.reshape(-1, 1)

r_list, v_list = ionsim.calculate_trajectory(
	r0, v0, 
	charge, mass, 
	20000, 0, 20, 
	force
)

assert(len(r_list) == 20000)
assert(len(v_list) == 20000)
assert(np.abs(np.array([r[0] for r in r_list]).max() - 1) < 0.005)
assert(np.abs(np.array([np.linalg.norm(v[0]) for v in v_list]).max() - 1) < 0.005)