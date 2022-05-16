import numpy as np

#values of exit and entrance pressure loss coefficient for turbulent flow, Re = 10 000
Kc = []
for sigma in range(0,10): #Kc seems linear
Kc =  0.51 - (0.41/1) * sigma
print(Kc)

Ke = [[0,1],
      [0.1,0.8],
      [0.2,0.62],
      [0.3,0.46],
      [0.4,0.32],
      [0.5,0.2],
      [0.6,0.1],
      [0.7,0.03],
      [0.8,-0.03],
      [0.9,-0.07],
      [1.-0.1]]

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