from models.freely_jointed_chain import FreelyJointedChain as FJC
from models.random_walks import RandomWalk as RW
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def calculate_end_to_end_distance(position):
    return np.linalg.norm(position[:, -1])

def calculate_statistical_obs(distances):
    statistical_obs = {}
    statistical_obs['rms'] = np.linalg.norm(distances)/len(distances)
    statistical_obs['variance'] = statistical_obs['rms']**2 - (np.sum(distances)/len(distances))**2
    statistical_obs['rmsf'] = statistical_obs['variance']*len(distances)/(len(distances)-1)
    statistical_obs['std_error_est'] = statistical_obs['variance']/(len(distances)-1)
    return statistical_obs


rw_sim = RW(nsteps=100, self_avoiding=True, can_walk_backwards=False)
fjc_sim = FJC(nsteps=100, self_avoiding=True, can_walk_backwards=False)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(rw_sim.position[0,:], rw_sim.position[1,:], rw_sim.position[2,:], label='Random walk')
ax.plot(fjc_sim.position[0,:], fjc_sim.position[1,:], fjc_sim.position[2,:], label='Freely jointed chain')
plt.legend()
plt.show()
