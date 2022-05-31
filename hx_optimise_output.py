import hx_functions as hxf
from hx_classes import water
import numpy as np


# hxf.brute_opt(n=1)
# h_w = water(55,45)
# c_w = water(20,25)

# hx_name = 'hx_y2018_p2022'

# hx = hxf.heat_exchangers(hx_name)[hx_name]

# hydraulic = hxf.hydraulic_design(0.5,0.5,h_w,c_w,hx,K_hot = 1.8,K_cold = 1)

# print(hydraulic['m_h'],hydraulic['m_c'])


#print(hxf.brute_opt_2())

k_array = [1,1,1,1,1,1,1,1,1,1,1]

print(k_array)

hxf.predict_hx(data = 'moodle',k_array=k_array)

hxf.predict_hx(data = '2022',k_array=k_array)
