from numba import jit
import numpy as np
import random


@jit(forceobj=True)
def check_if_crossing(current_step, self_avoiding_radius, position):
    for n in range(current_step+1):
        diff = position[current_step+1, :] - position[n, :]
        if self_avoiding_radius > np.sqrt(diff.dot(diff)):
            return True
    return False


class FreelyJointedChain:
    """Simulates a polymer using a freely jointed chain in 3D."""

    def __init__(self, nsteps, length=1, self_avoiding_radius=0.01,
                    self_avoiding=False, can_walk_backwards=True):

        self.nsteps = nsteps
        self._lenght = length
        self._self_avoiding = self_avoiding
        self._self_avoiding_radius = self_avoiding_radius
        self._can_walk_backwards = can_walk_backwards

        self.position = np.zeros((self.nsteps+1, 3), dtype='float64')

        self.successfully_self_avoidning = True
        self.step_number_when_breaking = 0

        self._run_simulation()

    def __str__(self):
        return 'FreelyJointedChain'

    def _run_simulation(self):

        for n in range(self.nsteps):
            new_vector = self._generate_new_vector(n)
            self.position[n+1, :] = self.position[n, :] + new_vector
            if self._self_avoiding:
                crosses_itself = self._check_if_crossing_itself(current_step=n)
                if crosses_itself:
                    self.successfully_self_avoidning = False
                    self.step_number_when_breaking = n
                    break

    def _generate_new_vector(self, n):
        """Generates new vector, if can_walk_backwards is False then checks that
        the angle between the new and former vector is not larger than
        arcsin(self.self_avoiding_radius/self.lenght)."""
        # Spherical coordinates
        phi = random.uniform(0, 2*np.pi)
        costheta = random.uniform(-1, 1)

        theta = np.arccos(costheta)

        x = self._lenght * np.sin(phi) * np.sin(theta)
        y = self._lenght * np.cos(phi) * np.sin(theta)
        z = self._lenght * np.cos(theta)

        new_vector = np.array([x, y, z])

        if not self._can_walk_backwards:
            # Uses a*b/|a||b|=cos(theta)
            a = new_vector
            b = self.position[n, :]
            angle_between_vectors = np.arccos(a.dot(b)/(np.linalg.norm(a)*np.linalg.norm(b)))
            if angle_between_vectors < np.arcsin(self._self_avoiding_radius/self._lenght):
                return self._generate_new_vector(n)
        return new_vector

    def _check_if_crossing_itself(self, current_step):
        """Returns True if walk crosses itself else False"""
        return check_if_crossing(current_step , self._self_avoiding_radius, self.position)
