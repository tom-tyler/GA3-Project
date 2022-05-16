import numpy as np
from hx_functions import KcKe

sigma = 0.05

print(KcKe(sigma))
#values of exit and entrance pressure loss coefficient for turbulent flow, Re = 10 000

mdot_dp_h = [[0.4583,   0.1333],
             [0.4236,   0.1756],
             [0.4010,   0.2024],
             [0.3611,   0.2577],
             [0.3125,   0.3171],
             [0.2639,   0.3633],
             [0.2222,   0.4233],
             [0.1597,   0.4784],
             [0.1181,   0.5330],
             [0.0694,   0.5715]] #[mass flow rate,pressure difference] for hot side

mdot_dp_c = [[0.5833,   0.1113],
             [0.5083,   0.2157],
             [0.4750,   0.2538],
             [0.4250,   0.3168],
             [0.3792,   0.3613],
             [0.3417,   0.4031],
             [0.2958,   0.4511],
             [0.2583,   0.4846],
             [0.2125,   0.5181],
             [0.1708,   0.5573]] #[mass flow rate, pressure difference] for cold side
