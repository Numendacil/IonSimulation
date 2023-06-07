import ionsim 
import numpy as np

from matplotlib import pyplot as plt

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

r = r0
v = v0

error = np.array([])

for i in range(0, 100, 10):
	r_list, v_list = ionsim.calculate_trajectory(
		r, v, 
		charge, mass, 
		1000, i, i + 10, 
		force
	)

	print(np.abs(np.array([np.linalg.norm(r[0]) for r in r_list]).max() - 1), i)

	r = r_list[-1]
	v = v_list[-1]

	error = np.append(error, np.abs([np.linalg.norm(r[0]) - 1 for r in r_list]))

plt.plot(np.linspace(0, 100, error.shape[0]), error)
plt.show()

# assert(len(r_list) == 20000)
# assert(len(v_list) == 20000)
# assert(np.abs(np.array([r[0] for r in r_list]).max() - 1) < 0.005)
# assert(np.abs(np.array([np.linalg.norm(v[0]) for v in v_list]).max() - 1) < 0.005)