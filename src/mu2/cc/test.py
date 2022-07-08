import sys
import numpy as np

from cscatter import k3cotdelta_py

sys.path.append('/Users/danielodell/6Li')
from system import System
import utility as util
import constants as const

R = 5.5

r, wr = util.log_mesh(0, 10*const.BETA6, 2000)

s = System(R, r, wr)

x = s.k3cotd_gen(s.ks, 0, 0)
xp = np.array([k3cotdelta_py(k, s.v_tilde, s.q, s.wq, 10*2/R, const.MASS) for k in s.ks])

print((x - xp) / x)