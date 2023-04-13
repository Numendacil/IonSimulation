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

# All vector inputs are given in shape of (3, N) rather than (N, 3) to improve speed. Use np.transpose() if necessary

# Acceleration of ions at given time with given location and velocity
def accl(r, v, t):
	return - 2.0 * r

N = 2
r0 = np.array([[0, 1, 0], [1, 0, 0]]).transpose()
v0 = np.random.rand(3, N)
step = 100
time_start = 0
time_end = 1

r_list, v_list = ionsim.calculate_trajectory(
	r0, v0,
	step,
	time_start,
	time_end,
	accl
)
```