from geneal.genetic_algorithms import ContinuousGenAlgSolver
from geneal.applications.fitness_functions.continuous import fitness_functions_continuous
from hx_classes import HX
import hx_functions as hxf

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

hx = HX(tube_number = 13,
        baffle_number = 14,
        pitch = 12e-3,
        tube_length = 362e-3,
        shell_length = 450e-3,
        baffle_gap = 14e-3,
        baffle_type = 'across_c',
        tube_layout='t',
        shell_passes=1,
        nozzle_bore=20e-3,
        crossflow_tube_fraction = 1,
        bypass_area_fraction = 0,
        seal_strips = 0,
        crossflow_rows = 4.5,
        tube_bundle_diameter= 56e-3)

def brute_opt(hx):
    #brute force optimisation with a few checks to eliminate cases as early as possible

    baffle_type = 'across_c'
    tube_layout = 't'
    number_iterations = 0

    for tube_number in range(1,15):
        pitch_max = round(hx.shell.d_inner/tube_number)
        for tube_length in range(0,350):
            if tube_number * tube_length <= 3500:
                for shell_length in range(0,350):
                    for shell_passes in range(0,4):
                        for baffle_number in range(0,60):
                            for baffle_gap in range(0,64):
                                for seal_strips in range(0,20):
                                    if hx.total_mass() <= 1.1:
                                        for pitch in range(10,pitch_max):
                                            print(number_iterations)
                                            crossflow_rows = hx.shell.d_inner/pitch
                                            hx = HX(tube_number = tube_number,
                                                    baffle_number = baffle_number,
                                                    pitch = pitch/1000,
                                                    tube_length = tube_length/1000,
                                                    shell_length = shell_length/1000,
                                                    baffle_gap = baffle_gap/1000,
                                                    baffle_type = baffle_type,
                                                    tube_layout = tube_layout,
                                                    shell_passes = shell_passes,
                                                    nozzle_bore=20e-3,
                                                    crossflow_tube_fraction = 1,
                                                    bypass_area_fraction = 0,
                                                    seal_strips = seal_strips,
                                                    crossflow_rows = crossflow_rows,
                                                    tube_bundle_diameter= 56e-3)
                                            number_iterations += 1


brute_opt(hx)
