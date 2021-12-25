import math

import numpy as np
import random as rnd

from matplotlib import pyplot as plt

r_list = [1]
random_parameters = [3, 4, 128]     # [a, c, m]


def generate_random_number(random_function_nr, range_randint):
    if random_function_nr == 1:
        return rnd.choice(range_randint)
    a = random_parameters[0]
    c = random_parameters[1]
    m = random_parameters[2]
    r = (a*r_list[-1] + c) % m
    r_list.append(r)
    return r//((m+3)//4)


def create_walk(number_of_steps, random_function_nr=1, can_cross_itself=True, can_walk_backwards=True):
    x_pos = [0]
    y_pos = [0]
    last_rand_int = 0
    default_range_randint = [1, 2, 3, 4]
    a={1:2,2:1,3:4,4:3}
    for n in range(number_of_steps-1):
        range_randint = default_range_randint.copy()
        if not can_walk_backwards:
            try:
                range_randint.remove(a[last_rand_int])
            except KeyError:
                pass
        rand_int = generate_random_number(random_function_nr, range_randint)
        last_rand_int = rand_int
        if rand_int == 1:
            x_pos.append(x_pos[-1]+1)
            y_pos.append(y_pos[-1])
        elif rand_int == 2:
            x_pos.append(x_pos[-1]-1)
            y_pos.append(y_pos[-1])
        elif rand_int == 3:
            y_pos.append(y_pos[-1]+1)
            x_pos.append(x_pos[-1])
        elif rand_int == 4:
            y_pos.append(y_pos[-1]-1)
            x_pos.append(x_pos[-1])

        if not can_cross_itself:
            if check_if_walk_crosses_itself(x_pos, y_pos):
                return x_pos, y_pos, 'crosses itself'

    return x_pos, y_pos, 'does not cross itself'


def plot_walk(*args, title='Random walk'):
    plt.clf()
    for i in range(0, len(args), 2):
        plt.scatter(args[i], args[i+1])
        plt.plot(args[i], args[i+1])
    plt.title(title)
    plt.show()


def calculate_distance(x_pos, y_pos):
    return math.sqrt(x_pos[-1]**2+y_pos[-1]**2)


def calculate_root_mean_squared(lst, n):
    return np.sqrt(np.sum(lst**2))/n


def calculate_root_mean_squared_fluctuation(lst, n):
    a = calculate_root_mean_squared(lst, n)**2
    b = calculate_root_mean_squared(lst**2, n)
    return np.sqrt(b-a)*n/(n-1)


def calculate_standard_error(lst, n):
    return [x/np.sqrt(n) for x in lst]


def check_if_walk_crosses_itself(x_pos, y_pos):
    """Returns True if walk crosses itself else False."""
    if (x_pos[-1], y_pos[-1]) in zip(x_pos[:-1], y_pos[:-1]):
        return True
    return False


def assignment_1a():
    for n in [10, 100, 1000]:
        x_pos, y_pos, a = create_walk(n)
        plot_walk(x_pos, y_pos, title='Random walk with '+str(n)+' steps')


def assignment_1b():
    # for a in range(1, 100, 10):
    #     random_parameters[0] = a
    #     x_pos, y_pos, z = create_walk(10, random_function_nr=2)
    #     plot_walk(x_pos, y_pos, title='Random walk with a = '+str(a))

    # for c in range(10):
    #     random_parameters[1] = c
    #     x_pos, y_pos, z = create_walk(10, random_function_nr=2)
    #     plot_walk(x_pos, y_pos, title='Random walk with c = '+str(c))

    # x_pos, y_pos = create_walk(100, random_function_nr=2)
    # plot_walk(x_pos, y_pos, title='Random walk with m = 128')
    # random_parameters[2]=101
    # x_pos, y_pos,z = create_walk(100, random_function_nr=2)
    # plot_walk(x_pos, y_pos, title='Random walk with m = 101')

    x_pos1, y_pos1, z = create_walk(15)
    random_parameters[2] = 101
    x_pos2, y_pos2, z = create_walk(100, random_function_nr=2)
    random_parameters[2] = 128
    x_pos3, y_pos3, z = create_walk(100, random_function_nr=2)
    random_parameters[0] = 10
    x_pos4, y_pos4, z = create_walk(15, random_function_nr=2)
    random_parameters[0] = 10
    x_pos5, y_pos5, z = create_walk(15, random_function_nr=2)
    plt.clf()
    plt.scatter(x_pos1, y_pos1)
    plt.plot(x_pos1, y_pos1, label='random.choice()')
    plt.scatter(x_pos2, y_pos2)
    plt.plot(x_pos2, y_pos2, label='other function with m=101')
    plt.title('Paths for different random number generators')
    plt.legend()
    plt.show()
    plt.clf()
    plt.scatter(x_pos1, y_pos1)
    plt.plot(x_pos1, y_pos1, label='random.choice()')
    plt.scatter(x_pos3, y_pos3)
    plt.plot(x_pos3, y_pos3, label='other function with m=128')
    plt.title('Paths for different random number generators')
    plt.legend()
    plt.show()

    plt.clf()
    plt.scatter(x_pos2, y_pos2)
    plt.plot(x_pos2, y_pos2, label='other function with m=101')
    plt.scatter(x_pos3, y_pos3)
    plt.plot(x_pos3, y_pos3, label='other function with m=128')
    plt.title('Paths for different random number generators')
    plt.legend()
    plt.show()


def assignment_1c():
    plt.figure()
    step_numbers = [x for x in range(20, 1000, 150)]
    number_of_walks = 1100
    x = np.arange(1, number_of_walks, 20)
    distance_v_step_numbers = {}
    for n in step_numbers:
        distance_v_step_numbers[n] = np.zeros(number_of_walks)
        for i in range(number_of_walks):
            x_pos, y_pos, z = create_walk(n)
            distance_v_step_numbers[n][i] = calculate_distance(x_pos, y_pos)
        y = np.array([np.sqrt(calculate_root_mean_squared(distance_v_step_numbers[n][:i]**2, i)) for i in x])
        plt.loglog(x, y, label="N = " + str(n))
        
    
    plt.title('RMS distance vs number of walks for different N')
    plt.legend()
    plt.xlabel("number of walks")
    plt.ylabel("rms distance")
    plt.show()

    plt.figure()
    for n in step_numbers:
        y = np.array([calculate_root_mean_squared_fluctuation(distance_v_step_numbers[n][:i], i) for i in x])
        plt.loglog(x, y, label="N = " + str(n))
    plt.title('RMS distance fluctuation vs number of walks for different N')
    plt.legend()
    plt.xlabel("number of walks")
    plt.ylabel("rmsf")
    plt.show()

    plt.figure()
    for n in step_numbers:
        y = np.array([calculate_root_mean_squared_fluctuation(distance_v_step_numbers[n][:i], i)/np.sqrt(i-1) for i in x])
        plt.loglog(x, y, label="N = " + str(n))
    plt.title('Standard error vs number of walks for different N')
    plt.legend()
    plt.xlabel("number of walks")
    plt.ylabel("standard error")
    plt.show()

    # Plot with errorbars
    #plt.clf()
    #plt.errorbar(step_numbers, root_mean_squared, yerr=standard_errors)
    #plt.title('RMS distance vs number of steps, with errorbars')
    #plt.show()


def assignment_1d():
    step_numbers = [a for a in range(1,40)]
    successful_walks = {}
    number_of_tries = 100

    # If can walk backwards
    for n in step_numbers:
        successful_walks[n] = 0
        for i in range(1, number_of_tries):
            x_pos, y_pos, outcome = create_walk(n, can_cross_itself=False)
            if outcome == 'does not cross itself':
                successful_walks[n] += 1
    successful_walk_fractions = [x/number_of_tries for x in successful_walks.values()]

    a =0.35

    plt.figure()
    plt.scatter(step_numbers, successful_walk_fractions, label='fraction of successful walks')
    plt.plot(step_numbers, np.exp([-a*x for x in step_numbers]), label='e^-'+str(a)+'x')
    plt.title('Fraction of successful walks vs number of steps')
    plt.legend()
    plt.show()

    # If can't walk backwards
    for n in step_numbers:
        successful_walks[n] = 0
        for i in range(number_of_tries):
            x_pos, y_pos, outcome = create_walk(n, can_cross_itself=False, can_walk_backwards=False)
            if outcome == 'does not cross itself':
                successful_walks[n] += 1
    successful_walk_fractions = [x / number_of_tries for x in successful_walks.values()]
    plt.figure()
    plt.scatter(step_numbers, successful_walk_fractions, label='fraction of successful walks')
    plt.plot(step_numbers, [1.07**(-x) for x in step_numbers], label='1.07^-x')
    plt.title('Fraction of successful walks vs number of steps')
    plt.legend()
    plt.show()


def assignment_1e():
    # number_of_steps = 500
    # number_of_walks = 100
    # x = np.arange(1, number_of_walks, 20)
    # scenarios = {1: ['non self avoiding', True, True], 2: ['self avoiding', False, True], 3: ['self avoding and can\'t walk backwards', False, False]}
    # plt.figure()
    # for m in scenarios:
    #     z = np.zeros(number_of_walks)
    #     for i in range(number_of_walks):
    #         x_pos, y_pos, outcome = create_walk(number_of_steps, can_cross_itself=scenarios[m][1], can_walk_backwards=scenarios[m][2])
    #         if outcome == 'does not cross itself':
    #            z[i] = calculate_distance(x_pos, y_pos)
    #     y = np.array([np.sqrt(calculate_root_mean_squared(z[:i]**2, i)) for i in x])
    #     plt.loglog(x, y, label=scenarios[m][0])
    # plt.legend()
    # plt.title('RMS distance vs number of walks averaged over, N = 500')
    # plt.xlabel('Number of walks')
    # plt.ylabel('RMS distance')
    # plt.show()

    number_of_steps = 50
    number_of_walks = 100
    scenarios = {1: ['non self avoiding', True, True], 2: ['self avoiding', False, True], 3: ['self avoding and can\'t walk backwards', False, False]}
    plt.figure()
    nn = np.arange(5, number_of_steps, 1)
    for m in scenarios:
        y = []
        for n in nn:     
            z = np.zeros(number_of_walks)
            for i in range(number_of_walks):
                x_pos, y_pos, outcome = create_walk(n, can_cross_itself=scenarios[m][1], can_walk_backwards=scenarios[m][2])
                if outcome == 'does not cross itself':
                    z[i] = calculate_distance(x_pos, y_pos)
            y.append(np.sqrt(calculate_root_mean_squared(z**2, len(z))))
        plt.scatter(nn, y, label=scenarios[m][0])
    plt.legend()
    plt.title('RMS distance vs number of steps')
    plt.xlabel('Number of steps')
    plt.ylabel('RMS distance')
    plt.show()

def assig_1_nice_plots():
    for n in [1,2]:
        x_pos, y_pos, a = create_walk(1000, random_function_nr=n)
        plot_walk(x_pos, y_pos, title='Random walk with 1000 steps for function number '+str(n))

# Function calls
# assignment_1a()
# assignment_1b()
# assignment_1c()
# assignment_1d()
assignment_1e()
# assig_1_nice_plots()