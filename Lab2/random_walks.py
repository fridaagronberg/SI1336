import math

import numpy as np
import random as rnd

from matplotlib import pyplot as plt

r_list = [1]
random_parameters = [3, 4, 128]     # [a, c, m]


def generate_random_number(random_function_nr):
    if random_function_nr == 1:
        return rnd.randint(1, 4)
    a = random_parameters[0]
    c = random_parameters[1]
    m = random_parameters[2]
    r = (a*r_list[-1] + c) % m
    r_list.append(r)
    return r//((m+3)//4)


def create_walk(number_of_steps, random_function_nr=1, can_cross_itself=True, can_walk_back=True):
    x_pos = np.zeros(number_of_steps, dtype=np.int64)
    y_pos = np.zeros(number_of_steps, dtype=np.int64)

    for n in range(number_of_steps-1):
        rand_int = generate_random_number(random_function_nr)
        if rand_int == 1:
            x_pos[n+1] = x_pos[n]+1
        elif rand_int == 2:
            x_pos[n+1] = x_pos[n]-1
        elif rand_int == 3:
            y_pos[n+1] = y_pos[n]+1
        elif rand_int == 4:
            y_pos[n+1] = y_pos[n]-1

        if not can_cross_itself:
            check_if_walk_crosses_itself(x_pos, y_pos)

    return x_pos, y_pos


def plot_walk(*args, title='Random walk'):
    plt.clf()
    for i in range(0, len(args), 2):
        plt.plot(args[i], args[i+1])
    plt.title(title)
    plt.show()


def calculate_distance(x_pos, y_pos):
    return math.sqrt(x_pos[-1]**2+y_pos[-1]**2)


def calculate_root_mean_squared(array, n):
    return np.sqrt(np.dot(array, array)/n)


def calculate_root_mean_squared_fluctuation(array, n):
    return np.sqrt((np.dot(array, array)+np.sum(array)**2)*n/(n-1))


def calculate_standard_error(array, n):
    return array/np.sqrt(n)


def check_if_walk_crosses_itself(x_pos, y_pos):
    pass


def assignment_1a():
    for n in [10, 100, 1000]:
        x_pos, y_pos = create_walk(n)
        plot_walk(x_pos, y_pos, title='Random walk with '+str(n)+' steps')


def assignment_1b():
    # for a in range(1, 100, 10):
    #     random_parameters[0] = a
    #     x_pos, y_pos = create_walk(10, random_function_nr=2)
    #     plot_walk(x_pos, y_pos, title='Random walk with a = '+str(a))

    # for c in range(10):
    #     random_parameters[1] = c
    #     x_pos, y_pos = create_walk(10, random_function_nr=2)
    #     plot_walk(x_pos, y_pos, title='Random walk with c = '+str(c))

    x_pos, y_pos = create_walk(100, random_function_nr=2)
    plot_walk(x_pos, y_pos, title='Random walk with m = 128')
    random_parameters[2]=101
    x_pos, y_pos = create_walk(100, random_function_nr=2)
    plot_walk(x_pos, y_pos, title='Random walk with m = 101')


def assignment_1c():
    step_numbers = [x for x in range(1, 1000, 50)]
    number_of_iterations = 100
    distance_v_step_numbers = {}
    root_mean_squared = []
    root_mean_squared_fluctuation = []
    for n in step_numbers:
        distance_v_step_numbers[n] = np.zeros(number_of_iterations)
        for i in range(number_of_iterations):
            x_pos, y_pos = create_walk(n)
            distance_v_step_numbers[n][i] = calculate_distance(x_pos, y_pos)
        root_mean_squared.append(calculate_root_mean_squared(distance_v_step_numbers[n], number_of_iterations))
        root_mean_squared_fluctuation.append(calculate_root_mean_squared_fluctuation(distance_v_step_numbers[n], number_of_iterations))
    standard_errors = calculate_standard_error(root_mean_squared_fluctuation, number_of_iterations)

    # plt.clf()
    # plt.plot(step_numbers, root_mean_squared)
    # plt.title('RMS distance vs number of steps, with '+str(number_of_iterations)+' iterations for each N')
    # plt.show()
    #
    # plt.clf()
    # plt.plot(step_numbers, root_mean_squared_fluctuation)
    # plt.title('RMS distance fluctuation vs number of steps, with ' + str(number_of_iterations) + ' iterations for each N')
    # plt.show()
    #
    # plt.clf()
    # plt.plot(step_numbers, standard_errors)
    # plt.title('Standard error vs number of steps, with ' + str(number_of_iterations) + ' iterations for each N')
    # plt.show()

    # Plot with errorbars
    # plt.clf()
    # plt.errorbar(step_numbers, root_mean_squared, yerr=standard_errors)
    # plt.title('RMS distance vs number of steps, with errorbars')
    # plt.show()


# Function calls
# assignment_1a()
# assignment_1b()
# assignment_1c()
