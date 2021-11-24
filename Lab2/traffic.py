import math
import random

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


class TrafficSimulation:
    def __init__(self, road_length=50, number_of_cars=25, v_max=2, p=0.5,
                 number_of_timesteps=100, initial_pos=None, initial_vel=None, animate=False):
        self.road_length = road_length
        self.number_of_cars = number_of_cars
        self.v_max = v_max
        self.p = p
        self.number_of_timesteps = number_of_timesteps
        self.animate = animate

        self.pos = np.zeros((self.number_of_timesteps, self.number_of_cars))
        self.vel = np.zeros((self.number_of_timesteps, self.number_of_cars))
        self.flow = 0
        self.equilibrium_stepnr = 0

        self._set_inital_conditions(initial_pos, initial_vel)
        self._run()

    def _set_inital_conditions(self, initial_pos, initial_vel):
        if initial_pos:
            self.pos[0,:] = initial_pos
        else:
            if self.road_length >= self.number_of_cars:
                self.pos[0, :] = np.arange(self.number_of_cars)
            else:
                raise ValueError('road length has to be greater or equal to number of cars')
        if initial_vel:
            self.vel[0] = initial_vel
        else:
            self.vel[0, :] = self.v_max - 1

    def _run(self):
        def _calculate_distance(x_i, x_j):  # j = i+1
            return (x_j-x_i+self.road_length)%self.road_length

        def _randomly_lower_vel():
            if random.random() <= self.p:
                return True
            return False

        for n in range(1, self.number_of_timesteps):
            for i in range(self.number_of_cars):
                d = _calculate_distance(self.pos[n - 1, i], self.pos[n - 1, (i+1)%self.number_of_cars])

                self.vel[n,i] = self.vel[n-1,i]
                if self.vel[n,i] < self.v_max:
                    self.vel[n,i] += 1
                if self.vel[n,i] >= d:
                    self.vel[n,i] = d-1
                if self.vel[n,i] > 0 and _randomly_lower_vel():
                    self.vel[n,i] -= 1
                self.pos[n,i] = self.pos[n-1,i]+self.vel[n,i]
            if self.animate:
                self._animate(self.pos[n, :])

    def show(self):
        plt.clf()
        for car_num in range(self.number_of_cars):
            plt.scatter(self.pos[:, car_num] % self.road_length, np.arange(self.number_of_timesteps), marker='.')
        plt.xlabel('Position')
        plt.ylabel('Time')
        plt.title('Position and time for '+str(self.number_of_cars)+' cars and '+str(self.number_of_timesteps)+' steps')
        plt.show()

    def _animate(self, pos):
        angles = (2*math.pi*pos % self.road_length)/self.road_length
        fig, ax = plt.subplots()
        xdata, ydata = [], []
        ln, = plt.plot([], [], 'ro')

        def init():
            ax.set_xlim(-1.2, 1.2)
            ax.set_ylim(-1.2, 1.2)
            return ln,

        def update(angles):
            xdata.append(np.cos(angles))
            ydata.append( np.sin(angles))
            ln.set_data(xdata, ydata)
            return ln,

        ani = FuncAnimation(fig, update(angles), frames=np.linspace(0, 2 * np.pi, 128),
                            init_func=init, blit=True)
        plt.show()


sim = TrafficSimulation(number_of_timesteps=50, number_of_cars=10, animate=True)


