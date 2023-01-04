from typing import Callable
from .utility import ft_matrix_gen, Mesh
from .counterterm import Counterterm

class Interaction:
    def __init__(self,
        long_range_potential: Callable[[float, float], float],
        xterm: Counterterm,
        rmesh: Mesh,
        scheme: str = 'nonlocal'
    ):
        self.counterterm = xterm
        self.R = self.counterterm.R # R
        self.scheme = scheme

        if scheme == 'local' or scheme == 'semilocal':
            self.long_rang_potential = lambda r: long_range_potential(r, self.R)
        elif scheme == 'nonlocal':
            assert self.counterterm.nonloc, '''Counterterm must be nonlocal for nonlocal-type Interaction.'''
            self.long_rang_potential = lambda r: long_range_potential(r, self.R/10)

        self.rmesh = rmesh


    def matrix_elements(self, qmesh, l, lp):
        v =  ft_matrix_gen(self.long_rang_potential, l, lp, qmesh.nodes,
            self.rmesh.nodes, self.rmesh.weights)
        if self.scheme == 'nonlocal':
            v *= self.counterterm.x_reg
        return v
        
