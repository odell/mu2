from typing import Callable
from .utility import ft_matrix_gen, Mesh
from .counterterm import Counterterm

class Interaction:
    def __init__(self,
        long_range_potential: Callable[[float, float], float],
        xterm: Counterterm,
        rmesh: Mesh
    ):
        self.counterterm = xterm
        self.R = self.counterterm.R # R
        if self.counterterm.nonloc:
            self.long_rang_potential = lambda r: long_range_potential(r, self.R/10)
        else:
            self.long_rang_potential = lambda r: long_range_potential(r, self.R)
        self.rmesh = rmesh


    def matrix_elements(self, qmesh, l, lp):
        v =   ft_matrix_gen(self.long_rang_potential, l, lp, qmesh.nodes,
            self.rmesh.nodes, self.rmesh.weights)
        if self.counterterm.nonloc:
            v *= self.counterterm.x_reg
        return v
        
