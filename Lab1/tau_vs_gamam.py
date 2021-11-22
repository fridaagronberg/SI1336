import math

from matplotlib import pyplot as plt

dt = 0.1
G = 9.8


class Observables:
    """ Class for storing observables for an oscillator """

    def __init__(self):
        self.time = []  # list to store time
        self.pos = []  # list to store positions
        self.vel = []  # list to store velocities
        self.energy = []  # list to store energy


class DampedHarmonic:
    def __init__(self, m=1, c=9, t0=0, theta0=1, dtheta0=0):
        self.m = m  # mass of the pendulum bob
        self.c = c  # c = g/L
        self.L = G / c  # string length
        self.t = t0  # the time
        self.theta = theta0  # the position/angle
        self.dtheta = dtheta0  # the velocity


def integrate_rk4(osc, obs, gamma):
    def accel(theta, dtheta, omega=3):
        return -omega ** 2 * theta - gamma * dtheta


    osc.t += dt
    a1 = accel(osc.theta, osc.dtheta) * dt
    b1 = osc.dtheta * dt
    a2 = accel(osc.theta + b1 / 2, osc.dtheta + a1 / 2) * dt
    b2 = (osc.dtheta + a1 / 2) * dt
    a3 = accel(osc.theta + b2 / 2, osc.dtheta + a2 / 2) * dt
    b3 = (osc.dtheta + a2 / 2) * dt
    a4 = accel(osc.theta + b3, osc.dtheta + a3) * dt
    b4 = (osc.dtheta + a3) * dt

    osc.dtheta += (a1 + 2 * a2 + 2 * a3 + a4) / 6
    osc.theta += (b1 + 2 * b2 + 2 * b3 + b4) / 6

    obs.pos.append(osc.theta)
    obs.vel.append(osc.dtheta)
    obs.time.append(osc.t)
    obs.energy.append(osc.m * osc.L * G * math.sin(osc.theta))


def run(tmax=30., title="tau vs gamma"):
    gammas=[x / 10.0 for x in range(1, 80)]
    taus = []
    for gamma in gammas:
        n = int(tmax / dt)
        osc = DampedHarmonic(theta0=math.pi * 0.5)
        obs = Observables()
        temp_tau = 0
        for it in range(n):
            integrate_rk4(osc, obs, gamma)
            if obs.pos[0]/math.e-0.1 < obs.pos[-1] < obs.pos[0]/math.e+0.1:
                temp_tau = obs.time[-1]
        taus.append(temp_tau)

        #plt.clf()
        #plt.title(title)
        #plt.plot(obs.time, obs.pos, label="Position")
        #plt.plot(obs.time, obs.vel, label="Velocity")
        #plt.plot(obs.time, obs.energy, label="Energy")
        #plt.xlabel('time')
        #plt.ylabel('observables')
        #plt.legend()
        #plt.show()
    y=[2/x for x in gammas]
    print(len(gammas), len(taus))
    plt.clf()
    plt.plot(gammas, taus, label='tau vs gamma')
    plt.plot(gammas,y, label='2/x')
    plt.xlabel('gamma')
    plt.ylabel('tau')
    plt.legend()
    plt.title('Relaxation time tau for different values on gamma')
    plt.show()
run()
