from numba import jit
import numpy as np

new_vectors = {
    1: np.array([1, 0, 0]),
    2: np.array([-1, 0, 0]),
    3: np.array([0, 1, 0]),
    4: np.array([0, -1, 0]),
    5: np.array([0, 0, 1]),
    6: np.array([0, 0, -1]),
    'last_rand_int': 0
}


@jit(nopython=True)
def check_if_crossing(current_step, position, new_pos):
    for n in range(current_step):
        if np.array_equal(position[n, :], new_pos):
            return True
    return False


class RandomWalk:
    """Simulates a polymer using a random walk in 3D"""
    def __init__(self, nsteps, length=1, self_avoiding=False, can_walk_backwards=True):
        self.nsteps = nsteps
        self._length = length
        self._self_avoiding = self_avoiding
        self._can_walk_backwards = can_walk_backwards

        self.position = np.zeros((self.nsteps + 1, 3), dtype='float64')

        self.successfully_self_avoidning = True
        self.step_number_when_breaking = self.nsteps

        self._run_simulation()

    def __str__(self):
        return 'RandomWalk'

    def _run_simulation(self):

        for n in range(self.nsteps):
            new_vector = self._generate_new_vector()
            new_pos = self.position[n, :] + new_vector
            if self._self_avoiding:
                crosses_itself = self._check_if_crossing_itself(n, new_pos)
                if crosses_itself:
                    self.successfully_self_avoidning = False
                    self.step_number_when_breaking = n
                    break
            self.position[n + 1, :] = new_pos


    def _generate_new_vector(self):
        random_number = np.random.randint(1, 7)
        if not self._can_walk_backwards:
            if new_vectors['last_rand_int'] == random_number:
                return self._generate_new_vector()
            new_vectors['last_rand_int'] = random_number
        return new_vectors[random_number]

    def _check_if_crossing_itself(self, current_step, new_pos):
        """Returns True if the random walk crosses itself, else False."""
        for i in range(0, current_step):
            if self.position[i, 0] == new_pos[0]:
                if self.position[i, 1] == new_pos[1]:
                    if self.position[i, 2] == new_pos[2]:
                        return True
        return False #check_if_crossing(current_step, self.position, new_pos)
