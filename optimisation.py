from geneal.genetic_algorithms import ContinuousGenAlgSolver
from geneal.applications.fitness_functions.continuous import fitness_functions_continuous

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

solver.solve()