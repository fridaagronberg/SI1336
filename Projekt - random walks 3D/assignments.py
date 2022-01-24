from models.freely_jointed_chain import FreelyJointedChain as FJC
from models.random_walks import RandomWalk as RW
import numpy as np
import matplotlib.pyplot as plt


def calculate_end_to_end_distance(position):
    return np.linalg.norm(position[:, -1])


def calculate_statistical_obs(distances):
    statistical_obs = {}
    statistical_obs['rms'] = np.linalg.norm(distances)/len(distances)
    statistical_obs['variance'] = statistical_obs['rms']**2 - (np.sum(distances)/len(distances))**2
    statistical_obs['rmsf'] = statistical_obs['variance']*len(distances)/(len(distances)-1)
    statistical_obs['std_error_est'] = statistical_obs['variance']/(len(distances)-1)
    return statistical_obs


def check_dependency_of_radius():
    radiuses = np.arange(0.01, 1, 0.01)# np.arange(0.001, 1, 0.01)
    res = {}
    nsteps = np.arange(5,100,30)  # np.arange(10000, 30000, 5000) #stepsize 500 sen
    number_of_sims = 1
    for nstep in nsteps:
        res[nstep] = []
        for r in radiuses:
            distances = []
            for i in range(number_of_sims):
                sim = FJC(self_avoiding_radius=r, nsteps=nstep, self_avoiding=True)
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.plot(sim.position[:, 0], sim.position[:, 1], sim.position[:, 2], label=r)
                plt.legend()
                plt.show()
                distances.append(calculate_end_to_end_distance(sim.position))
            res[nstep].append(calculate_statistical_obs(distances)['rms'])

    for nstep in nsteps:
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(radiuses, res[nstep], label='N = '+str(nstep))
        plt.legend()
        plt.xlabel('Self-avoiding radius')
        plt.ylabel('Distances')
        plt.show()

check_dependency_of_radius()

"""rw_sim = RW(nsteps=100, self_avoiding=False, can_walk_backwards=False)
fjc_sim = FJC(nsteps=100, self_avoiding=False, can_walk_backwards=False)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(rw_sim.position[0,:], rw_sim.position[1,:], rw_sim.position[2,:], label='Random walk')
ax.plot(fjc_sim.position[0,:], fjc_sim.position[1,:], fjc_sim.position[2,:], label='Freely jointed chain')
plt.legend()
plt.show()"""

