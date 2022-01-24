import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


class SimpleRelaxation:
    def __init__(self, expected_grid, length=10, steplength=1, left_boundary=10, right_boundary=10,
                 top_boundary=10, bottom_boundary=10, initial_grid=None):
        self.grid = np.zeros([length//steplength+1, length//steplength+1])
        self.length = length
        self.steplength = steplength
        self.expected_grid = expected_grid
        self.number_of_iterations = 0

        self._set_boundary(left_boundary, right_boundary, top_boundary, bottom_boundary)
        self._set_initial_values(initial_grid)

        # Sjukt obehagligt med det här är i init
        fig = plt.figure()
        ax = fig.add_subplot(111)
        self.im = ax.imshow(self.grid, cmap=None, interpolation='nearest')
        fig.colorbar(self.im)

        self.ani = animation.FuncAnimation(fig, self.calculate_potential, frames=500,
                                            interval=200, blit=True, repeat = False)
        #plt.show()

    def _set_boundary(self, left_boundary, right_boundary, top_boundary, bottom_boundary):
        self.grid[:, 0] = left_boundary
        self.grid[:, -1] = right_boundary
        self.grid[0, :] = top_boundary
        self.grid[-1, :] = bottom_boundary

    def _set_initial_values(self, inital_grid):
        if inital_grid.any() == None:
            pass
        else:
            self.grid[1:-1, 1:-1] = inital_grid

    def _simple_relaxtion_iteration(self):
        new_grid = np.zeros([self.length//self.steplength-2, self.length//self.steplength-2])
        for i in range(1, self.length//self.steplength):
            for j in range(1, self.length // self.steplength):
                new_grid[i-1, j-1] = 0.25*(self.grid[i+self.steplength, j] + self.grid[i-self.steplength, j]
                                            + self.grid[i, j-self.steplength] + self.grid[i, j+self.steplength])
        x = (new_grid - self.expected_grid) / self.expected_grid
        error = abs(max(x.min(), x.max(), key=abs))
        self.grid[1:-1, 1:-1] = new_grid
        return error

    def calculate_potential(self, framedata):
        error = self._simple_relaxtion_iteration()
        self.number_of_iterations += 1
        if error < 0.01:
            self.ani.save()
            self.ani.event_source.stop()
            print('Number of iterations required ', self.number_of_iterations)
        self.im.set_array(self.grid)
        return self.im,

initial_grid = 9*np.ones([9, 9])
SimpleRelaxation(expected_grid=initial_grid+1, initial_grid=initial_grid)
