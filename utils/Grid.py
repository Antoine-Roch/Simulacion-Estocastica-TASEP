import numpy as np
import copy
import time
from .InitialLey import *
from .MoveLey import *
from .LeyNumberParticules import *

class Grid:
    def __init__(self, N, initial_dist: Initial_Ley, n_step, move_ley: Move_Ley, ley_N: Ley_Number_Particules, 
                 x_lines = [0], y_lines = [], simul_type = "sequential", verif_mode = False):
        '''
        N: number of initial particules
        initial_dist: initial distribution of the particules along all the axis
        n_step: number of steps for the simulation
        x_lines: coordinates of the x_axis (list)
        y_lines: coordinates of the y_axis (list)
        move_type: standard (only right), left: right+left
        simul_type: sequential, parallel
        ley_N: ley that generates N, fixed, other (UPDATE)
        move_ley: ley of one particule
        '''

        assert len(x_lines) > 0, "No x axis"

        self.N_particules = N
        self.x_coordinates = np.array(x_lines)
        self.y_coordinates = np.array(y_lines)

        self.simul_type = simul_type
        self.ley_N = ley_N
        self.n_step = n_step
        self.move_ley = move_ley
        if len(y_lines) == 0: 
            self.check_good_move_ley = [simul_type, "uni"]
        elif len(y_lines) == 1 and len(x_lines) == 1:
            self.check_good_move_ley = [simul_type, "multi"]
        else:
            self.check_good_move_ley = [simul_type, "multi_HD"]

        self.positions = initial_dist(N, self.x_coordinates, self.y_coordinates)
        self.func_step = self._choose_step()
        self.positions_record = []
        self.verif_mode = verif_mode

        self._update_neighbours()

    def _step_sequential_bis(self):
        n_part = self.positions.shape[0]
        for k in range(n_part):
            self.positions[k] = np.copy(self.move_ley(self.positions[k], k, self.positions))

    def _update_neighbours(self):
        if self.check_good_move_ley[1] == "uni":
            self.dist = self.positions[1:] - self.positions[:-1]
        elif self.check_good_move_ley[1] == "multi":
            self.dist = []
            self.dist.append(self.positions[0][1:] - self.positions[0][:-1])
            self.dist.append(self.positions[1][1:] - self.positions[1][:-1])
            #WARNING: if two parallel axis are separated of 1, may not work
        elif self.check_good_move_ley[1] == "multi_HD":
            self.dist = [[], []]
            for pos in self.positions[0]:
                self.dist[0].append(pos[1:] - pos[:-1])
            for pos in self.positions[1]:
                self.dist[1].append(pos[1:] - pos[:-1])
        

    def _step_sequential(self):
        self.positions, self.dist = self.move_ley(self.positions, self.dist, x_axis = self.x_coordinates, y_axis = self.y_coordinates, n_axis = self.x_coordinates.shape[0] + self.y_coordinates.shape[0])

    def _step_parallel(self):
        pass

    def _choose_step(self):
        if self.simul_type == "sequential":
            return self._step_sequential
        elif self.simul_type == "parallel":
            return self._step_parallel

    def run_simulation(self):
        t = time.time()
        for k in range(self.n_step):
            self.func_step()
            self.N_particules, self.positions, self.dist = self.ley_N(self.positions, self.dist, self.check_good_move_ley[1], self.N_particules)
            self.positions_record.append(copy.deepcopy(self.positions))

            if self.verif_mode:
                
                if self.check_good_move_ley[1] == "multi":
                    A = (self.dist[0] == 0).any()
                    B = (self.dist[1] == 0).any()
                    C = np.sum(self.positions[0] == 0) + np.sum(self.positions[1] == 0) > 1

                    if A or B or C:
                        print(self.positions_record[-2])
                        print(self.positions_record[-1])
                        print(k, A, B, C)
                        break

                elif self.check_good_move_ley[1] == "uni":
                    A = (self.dist == 0).any()
                    if A:
                        print(self.positions_record[-2])
                        print(self.positions_record[-1])
                        print(k, A)
                        break
                
                elif self.check_good_move_ley[1] == "multi_HD":
                    for d in self.dist[0]:
                        if (d==0).any():
                            print(k, d)
                            break
                    for d in self.dist[1]:
                        if (d==0).any():
                            print(k, d)
                            break
            
        print(time.time() - t)

    def make_movie(self, save_path):
        pass