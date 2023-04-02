import pysim 
import numpy as np

r0 = np.array([[5, 0, 0]]).reshape(3, 1)
v0 = np.array([[0, 5, 0]]).reshape(3, 1)
def accl(r: np.ndarray, v: np.ndarray, t: float):
    return - r * 5.00

r_list, v_list = pysim.runge_kutta(r0, v0, 20000, 0, 20, accl)

with open("r_list.txt", "w") as file:
    for arr in r_list:
        file.write(np.array2string(arr))