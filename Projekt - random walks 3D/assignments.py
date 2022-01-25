from models.freely_jointed_chain import FreelyJointedChain as FJC
from models.random_walks import RandomWalk as RW
import numpy as np
import matplotlib.pyplot as plt


def calculate_end_to_end_distance(position):
    return np.linalg.norm(position[-1, :])


def calculate_statistical_obs(distances):
    statistical_obs = {}
    statistical_obs['rms'] = np.linalg.norm(distances)/np.sqrt(len(distances))
    statistical_obs['variance'] = statistical_obs['rms']**2 - (np.sum(distances)/len(distances))**2
    statistical_obs['rmsf'] = statistical_obs['variance']*len(distances)/(len(distances)-1)
    statistical_obs['std_error_est'] = statistical_obs['variance']/(len(distances)-1)
    return statistical_obs


def check_dependency_of_radius():
    radiuses = np.arange(0.01, 0.5, 0.005)# np.arange(0.001, 1, 0.01)
    res = {}
    nsteps = np.arange(1000, 6000, 1000) #stepsize 500 sen
    number_of_sims = 10
    for nstep in nsteps:
        res[nstep] = []
        for r in radiuses:
            distances = []
            for i in range(number_of_sims):
                sim = FJC(self_avoiding_radius=r, nsteps=nstep, self_avoiding=True)
                distances.append(calculate_end_to_end_distance(sim.position))
            res[nstep].append(calculate_statistical_obs(distances)['rms'])

    fig = plt.figure()
    ax = fig.add_subplot()
    for nstep in nsteps:
        ax.plot(radiuses, res[nstep], label='N = '+str(nstep))
    plt.legend()
    plt.xlabel('Self-avoiding radius')
    plt.ylabel('Distances')
    plt.show()

def check_dependency_of_nsteps():
    nsteps = np.arange(10, 50, 5)
    number_of_sims = 10
    for model, name in zip([FJC, RW], ['FreelyJointedChain', 'RandomWalk']):
        res = {'rms':[],'rmsf':[],'std_error_est':[]}
        res2 = {'rms':[],'rmsf':[],'std_error_est':[], 'numb_succesful_sims':[]}
        for nstep in nsteps:
            distances_all = []
            distances_success = []
            success_sims_counter = 0
            for i in range(number_of_sims):
                sim = model(nsteps=nstep, self_avoiding=True, can_walk_backwards=False)
                if sim.successfully_self_avoidning:
                    success_sims_counter += 1
                    distances_success.append(calculate_end_to_end_distance(sim.position))
                distances_all.append(calculate_end_to_end_distance(sim.position))
            stat_obs_all = calculate_statistical_obs(distances_all)
            res['rms'].append(stat_obs_all['rms'])
            res['rmsf'].append(stat_obs_all['rmsf'])
            res['std_error_est'].append(stat_obs_all['std_error_est'])
            stat_obs_success = calculate_statistical_obs(distances_success)
            res2['rms'].append(stat_obs_success['rms'])
            res2['rmsf'].append(stat_obs_success['rmsf'])
            res2['std_error_est'].append(stat_obs_success['std_error_est'])
            res2['numb_succesful_sims'].append(success_sims_counter)
        print(res)

        for i, key in enumerate(res.keys()):
            plt.figure(i+1)
            plt.plot(nsteps, res[key], label=name)
            plt.figure(i+4)
            plt.plot(nsteps, res2[key], label=name)
        plt.figure(7)
        plt.scatter(nsteps, res2['numb_succesful_sims'], label=name)


    for i, key in enumerate(['RMS', 'RMSF', 'Standard error estimate']):
        plt.figure(i+1)
        plt.ylabel(key)
        plt.legend()
        plt.xlabel('Number of steps')
        plt.figure(i+4)
        plt.ylabel(key)
        plt.legend()
        plt.xlabel('Number of steps')
    plt.figure(7)
    plt.ylabel('Number of succesful simulations')
    plt.legend()
    plt.xlabel('Number of steps')
    plt.show()

def check_dependency_of_step_length():
    plt.figure()
    step_lengths = np.arange(0.5, 3, 0.5)
    number_of_sims = 10
    for model, name in zip([FJC, RW], ['FreelyJointedChain', 'RandomWalk']):
        res = []
        for step_length in step_lengths:
            distances = []
            for i in range(number_of_sims):
                sim = model(nsteps=100, length=step_length, self_avoiding=True)
                distances.append(calculate_end_to_end_distance(sim.position))
            stat_obs = calculate_statistical_obs(distances)
            res.append(stat_obs['rms'])
        plt.plot(step_lengths, res, label=name)
    plt.legend()
    plt.xlabel('Step length')
    plt.ylabel('RMS distance')
    plt.show()

#check_dependency_of_radius()
#check_dependency_of_nsteps()
#check_dependency_of_step_length()


fjc_sim0 = FJC(nsteps=100, self_avoiding=False, can_walk_backwards=False)
fjc_sim1 = FJC(nsteps=500, self_avoiding=False, can_walk_backwards=False)
fjc_sim2 = FJC(nsteps=100, self_avoiding=True, can_walk_backwards=False)
fjc_sim3 = FJC(nsteps=500, self_avoiding=True, can_walk_backwards=False)
for i, sim in enumerate([fjc_sim0,fjc_sim1,fjc_sim2,fjc_sim3]):
    fig = plt.figure(i)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(sim.position[:,0], sim.position[:,1], sim.position[:,2])

plt.show()
