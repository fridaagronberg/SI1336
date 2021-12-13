import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


class LaplaceApp:
    def __init__(self, length=10, steplength=1, left_boundary=10, right_boundary=10,
                 top_boundary=10, bottom_boundary=10, inital_grid=None, method=_simple_relaxtion_iteration):
        self.grid = np.zeros(length//steplength, length//steplength)
        self.length = length
        self.steplength = steplength
        self._set_boundary(left_boundary, right_boundary, top_boundary, bottom_boundary)
        self._set_initial_values(inital_grid)
        self.method = method

    def _set_boundary(self, left_boundary, right_boundary, top_boundary, bottom_boundary):
        self.grid[:, 0] = left_boundary
        self.grid[:, -1] = right_boundary
        self.grid[0, :] = top_boundary
        self.grid[-1, :] = bottom_boundary

    def _set_initial_values(self, inital_grid):
        if inital_grid:
            self.grid[1:-1, 1:-1] = inital_grid
        else:
            pass

    def _simple_relaxtion_iteration(self):
        new_grid = np.zeros(self.length//self.steplength-1, self.length//self.steplength-1)
        for i in range(1, self.length//self.steplength):
            for j in range(1, self.length // self.steplength):
                new_grid[i-1, j-1] = 0.25*(self.grid[i+self.steplength, j] + self.grid[i-self.steplength, j]
                                            + self.grid[i, j-self.steplength] + self.grid[i, j+self.steplength])
        error = max(abs((new_grid - self.grid[1:-1, 1:-1]) / self.grid[1:-1, 1:-1]))
        self.grid[1:-1, 1:-1] = new_grid
        return error

    def calculate_potential(self, max_error=0.01):
        # Initialize animation?
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.imshow(v, cmap=None, interpolation='nearest')
        fig.colorbar(im)
        number_of_iterations = 0
        while True:
            error = self.method()
            number_of_iterations += 1
            if error < max_error:
                break
        print('Number of iterations required ', number_of_iterations)

