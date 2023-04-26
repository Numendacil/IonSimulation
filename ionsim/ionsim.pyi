from typing import Callable, List, Tuple

import numpy

def calculate_trajectory(
	init_r: numpy.ndarray, 
    init_v: numpy.ndarray, 
    charge: numpy.ndarray, 
    mass: numpy.ndarray, 
    step: int, 
    time_start: float, 
    time_end: float, 
    force: Callable[[numpy.ndarray, numpy.ndarray, float], numpy.ndarray]
) -> Tuple[List[numpy.ndarray],List[numpy.ndarray]]: ...
