import numpy as np

#values of exit and entrance pressure loss coefficient for turbulent flow, Re = 10 000
Kc = []
for sigma in range(0,10): #Kc seems linear
Kc =  0.51 - (0.41/1) * sigma])
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

mdot_dp = [[]]