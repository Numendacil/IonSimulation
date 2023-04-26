# IonSimulation

Currently WIP

### Install
Clone this repo and run
```sh
pip install ./IonSimulation
```

### Example usage

```python
import ionsim
import numpy as np

# External force on ions at given time with given location and velocity
def force(r, v, t):
	return - 2.0 * r

N = 2
r0 = np.array([[0, 1, 0], [1, 0, 0]])
v0 = np.random.rand(3, N)
charge = np.array([0])
mass = np.array([1])
step = 100
time_start = 0
time_end = 1

r_list, v_list = ionsim.calculate_trajectory(
	r0, v0,
	charge, mass,
	step,
	time_start,
	time_end,
	force
)
```