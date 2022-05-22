from geneal.genetic_algorithms import ContinuousGenAlgSolver
from geneal.applications.fitness_functions.continuous import fitness_functions_continuous
from hx_classes import HX
import hx_functions as hxf
import hx_design as hxd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

solver = ContinuousGenAlgSolver(n_genes=4, # number of variables defining the problem
                                fitness_function=fitness_functions_continuous(3), # fitness function to be maximized
                                pop_size=10,
                                max_gen=200,
                                mutation_rate=0.1,
                                selection_rate=0.6,
                                selection_strategy="roulette_wheel",
                                problem_type=float, # Defines the possible values as float numbers
                                variables_limits=(-10, 10) # Defines the limits of all variables between -10 and 10. 
                                                           # Alternatively one can pass an array of tuples defining the limits
                                                           # for each variable: [(-10, 10), (0, 5), (0, 5), (-20, 20)]
                                )

#solver.solve()




            
def optimiser(hx):
    best_array = []
    step_array = []
    tube_length_array = []
    step = 50

    optimised = brute_opt(hx,step)
    best_design = optimised[0]
    q_best = best_design[0]
    tube_length = optimised[1]
    best_array.append(q_best[0])
    tube_length_array.append(tube_length*20)
    step_array.append(step)

    #plt.plot(step_array, best_array)
    #plt.plot(step_array, tube_length_array)
    plt.show()

    return optimised


print(optimiser(HX))



#brute_opt(HX)

