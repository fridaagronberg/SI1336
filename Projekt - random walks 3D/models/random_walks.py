import numpy as np

new_vectors = {
1: np.array([1, 0, 0]),
2: np.array([-1, 0, 0]),
3: np.array([0, 1, 0]),
4: np.array([0, -1, 0]),
5: np.array([0, 0, 1]),
6: np.array([0, 0, -1]),
last_rand_int: 0
}

class RandomWalk:
    """Simulates a polymer using a random walk in 3D"""
    def __init__(self, nsteps, length=1, self_avoiding=False, can_walk_backwards=True):
        self.nsteps = nsteps
        self._length = length
        self._self_avoiding = self_avoiding
        self._can_walk_backwards = can_walk_backwards

        self.position = np.zeros((3, self.nsteps+1))

        self.successfully_self_avoidning = True

        self._run_simulation()

    def _run_simulation(self):

        for n in range(nsteps):
            new_vector = self._generate_new_vector(n):
            self.position[:, n+1] = self.position[:, n] + new_vector
            if self._self_avoiding:
                crosses_itself = self._check_if_crossing_itself(current_step=n)
                if crosses_itself:
                    self.successfully_self_avoidning = False
                    self.step_number_when_breaking = n
                    break

    def _generate_new_vector(self):
        random_number = np.randint(1, 7)
        if not self._can_walk_backwards:
            if new_vectors[last_rand_int] == random_number:
                return self._generate_new_vector()
            new_vectors[last_rand_int] = random_number
        return new_vectors[random_number]

    def _check_if_crossing_itself(current_step):
        """Returns True if the random walk crosses itself, else False."""
        for n in range(current_step):
            if self.position[:, n] == self.position[:, current_step]:
                return True
        return False
