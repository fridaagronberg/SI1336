import math
import random

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from numpy.lib.function_base import average


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
        self.flow = np.zeros(self.number_of_timesteps)
        self.average_flow = 0

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

            self.flow[n] = self.vel[n,:].sum()/self.road_length

            if self.animate:
                self._animate(self.pos[n, :])
        avg_flow = 0
        for n in range(20, self.number_of_timesteps):
            avg_flow += self.flow[n]
        self.average_flow = avg_flow/(self.number_of_timesteps-20)

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


# 2.2a)
# simulations = {}
# number_of_simulations = 100
# avg_flow = {}
# for num_cars in range(2,52,5):
#     simulations[num_cars] = []
#     avg_flow[num_cars] = 0
#     for i in range(number_of_simulations):
#         sim = TrafficSimulation(number_of_cars=num_cars)
#         simulations[num_cars].append(sim)
#         avg_flow[num_cars] += sim.average_flow
#     avg_flow[num_cars] = avg_flow[num_cars]/number_of_simulations
#
#
# plt.clf()
# plt.plot(np.arange(2,52,5)/50, avg_flow.values(), '.-')
# plt.xlabel('Density cars/road')
# plt.ylabel('Average flow')
# plt.title('Average flow vs. road density')
# plt.show()

def calculate_standard_error_estimate(lst):
    std = np.std(lst)
    return std/len(lst)

#2.b
#avg_flows = []
#errors = {}
#x=np.arange(2,20)
#for n in x:
#    sim = TrafficSimulation()
 #   avg_flows.append(sim.average_flow)
#
#for n in x:
#    errors[n] = calculate_standard_error_estimate(avg_flows[:n])
#
#plt.clf()
#plt.plot(errors.keys(), errors.values(), '.-')
#plt.xlabel('Number of simulations')
#plt.ylabel('Standard error estimate of flow')
#plt.title('Standard error estimate of flow vs number of simulations')
#plt.show()

# simulations = {}
# number_of_simulations = 100
# avg_flow = {}
# for length in range(25,100, 10):
#     simulations[length] = []
#     avg_flow[length] = 0
#     for i in range(number_of_simulations):
#         sim = TrafficSimulation(road_length=length)
#         simulations[length].append(sim)
#         avg_flow[length] += sim.average_flow
#     avg_flow[length] = avg_flow[length]/number_of_simulations
#
#
# plt.clf()
# plt.plot(25/np.arange(25,100, 10), avg_flow.values(), '.-')
# plt.xlabel('Density cars/road')
# plt.ylabel('Average flow')
# plt.title('Average flow vs. road density')
# plt.show()

#2.c
# n=10
# avg_flows=[]
# lengths = np.arange(25,1000, 50)
# for length in lengths:
#     avg_flow = 0
#     for i in range(n):
#         sim = TrafficSimulation(road_length=length, number_of_cars=int(0.25*length))
#         avg_flow += sim.average_flow
#     avg_flows.append(avg_flow/n)
# 
# plt.clf()
# plt.plot(lengths, avg_flows, '.-')
# plt.title("Average flow vs road length, for a given density of 0.25")
# plt.xlabel("Road length")
# plt.ylabel("Average flow")
# plt.show()

#2.d
max_vels = [1,2,5]
lengths = np.arange(25,200, 5)
avg_flows = {}
for v in max_vels:
    avg_flows[v] = []
    for lenght in lengths:
        sim = TrafficSimulation(road_length=lenght, v_max=v)
        avg_flows[v].append(sim.average_flow)

plt.clf()
plt.plot(25/lengths, avg_flows[1], '.-', label="v_max = 1")
plt.plot(25/lengths, avg_flows[2], '.-', label="v_max = 2")
plt.plot(25/lengths, avg_flows[5], '.-', label="v_max = 5")
plt.ylabel("Average flow")
plt.xlabel("Road density")
plt.title('Average flow vs. road density')
plt.legend()
plt.show()