import multiprocessing as mp

import matplotlib.pyplot as plt
import numpy as np

from utils import *


class ProcessPlotter:
	def __init__(self, queue_in: mp.Queue, queue_out: mp.Queue, frame_init : Frame, interval: float):
		self.queue_in = queue_in
		self.queue_out = queue_out
		self.fig, self.ax = plt.subplots(1, 2, figsize=(15, 10))

		for ax in self.ax:
			ax.set_xlim(-10.0, 10.0)
			ax.set_ylim(-10.0, 10.0)
			ax.set_aspect('equal')

		self.artists = (
			self.ax[0].plot(frame_init.r[:, 0], frame_init.r[:, 1], 'ro', animated=True)[0],
			self.ax[1].plot(frame_init.r[:, 0], frame_init.r[:, 2], 'ro', animated=True)[0],
			# self.ax[0].plot(frame_init.r[-1, 0], frame_init.r[-1, 1], 'bo', animated=True)[0],
			# self.ax[1].plot(frame_init.r[-1, 0], frame_init.r[-1, 2], 'bo', animated=True)[0],
		)
		self.bm = BlitManager(self.fig.canvas, self.artists)

		self.interval = interval
		self.count = 0

		plt.show(block=False)
		plt.pause(.1)

	def run(self):
		self.count += 1
		if self.queue_in.empty():
			return True
	
		f: Frame = self.queue_in.get()
		if f.timestamp < 1e-5:
			self.count = 0
		else:
			print(self.count * self.interval, f.timestamp)
		
			
		self.artists[0].set_data(f.r[:, 0], f.r[:, 1])
		self.artists[1].set_data(f.r[:, 0], f.r[:, 2])
		# self.artists[2].set_data(f.r[-1, 0], f.r[-1, 1])
		# self.artists[3].set_data(f.r[-1, 0], f.r[-1, 2])
		
		self.bm.update()

		return True

	def start(self,):
		print('starting plotter...')

		# timer = self.fig.canvas.new_timer(interval=int(self.interval * 1000))
		# timer.add_callback(self.run)
		# timer.start()

		while True:
			self.run()
			plt.pause(self.interval)


N = 50
# rand = [np.array([]), -100]
charge = np.ones(N)
mass = np.ones(N)
# mass[-1] = 0.8

if __name__ == "__main__":
	def force(r: np.ndarray, v: np.ndarray, t: float):

		aa = np.array([-0.00005, 0.0001, 0.4])
		q = np.array([0.1, 0, -0.1])

		gamma = (0.7 - t / 30 * 0.6) if t < 30 else 0.05
		# gamma = 0.5

		f: np.ndarray = - 2000 * r * (1 + 0.001 * np.sum(r * r, axis=0)) * charge.reshape(-1, 1) * (aa - q * np.cos(120 * t)) - gamma * v
		return f
	
		# E = 10
		# inter = 50
		# if 50 + inter >= t > 50:
		# 	f[:, 0] += E
		# elif 50 + 2 * inter >= t > 50 + inter:
		# 	f[:, 0] += -E
		# elif 50 + 3 * inter >= t > 50 + 2 * inter:
		# 	f[:, 0] += E
		# elif 50 + 4 * inter >= t > 50 + 3 * inter:
		# 	f[:, 0] += -E

		# global rand
		# if t - rand[1] > 0.05:
		# 	rand[0] = np.random.rand(N, 3) * 0.6 - 0.3
		# 	rand[1] = t
		# return f + rand[0]
	

	p = Producer(step=120, interval=0.04, batch=50)

	r0 = np.random.rand(N, 3) * 10 - 5
	# r0[:, -1] = np.array([30, 0, 0])
	v0 = np.zeros((N, 3))

	q1 = mp.Queue()
	q2 = mp.Queue()

	m = Message()
	m.command = CommandType.START
	m.r = r0
	m.v = v0
	m.mass = mass
	m.charge = charge
	m.force = force

	q1.put(m)
	q2.put(Frame(r0, v0, 0))
	plot = ProcessPlotter(q2, q1, Frame(r0, v0, 0), interval=0.04)
	proc = mp.Process(target=p.run, args=(q2, q1,), daemon=True)
	proc.start()

	plot.start()
