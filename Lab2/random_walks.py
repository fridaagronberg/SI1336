import numpy as np
import random as rnd

from matplotlib import pyplot as plt

r_list = [1]
random_parameters = [3, 4, 128]     # [a, c, m]


def generate_random_nummer(random_function_nr):
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
        rand_int = generate_random_nummer(random_function_nr)
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


def calculate_root_mean_squared():
    pass


def calculate_root_mean_squared_fluctuation():
    pass


def calculate_standard_error():
    pass


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

# Function calls
# assignment_1a()
assignment_1b()