"""

"""

import os
import numpy as np

# Defaults for when this file is run
TRACK_FILEPATH_DEFAULT = "input/debug_track_%03d.txt"
SIM_OUTPUT_FILEPATH_BASE_DEFAULT = "output/debug_nav_%03d.txt"
ERROR_OUTPUT_FILEPATH_BASE_DEFAULT = "output_evaluation/debug_nav_errors_%03d.txt"
MEAN_ERROR_OUTPUT_FILEPATH_DEFAULT = "output_evaluation/debug_nav_mean_errors.txt"
NUM_SIMULATIONS_DEFAULT = 1


# Creates files for the errors of simulations, and the mean errors of repeated simulations
def create_sim_error_files(track_filepath, sim_output_filepath_base, error_output_filepath_base, mean_error_output_filepath, num_simulations):

    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    # Get and save all the errors of each simulation to the track file
    all_errors = []
    for i in range(1, num_simulations+1):
        errors = []
        to_skip = False
        with open(track_filepath % i, 'r') as track_f:
            with open(sim_output_filepath_base % i, 'r') as sim_output_f:

                timesteps = int(track_f.readline())
                dimension = int(track_f.readline())
                init_state = np.array([float(x) for x in track_f.readline().split()])
                init_cov = np.array([[float(x) for x in track_f.readline().split()] for _ in range(dimension)])

                for t in range(timesteps):
                    
                    # Read ground truth from track file and estimate from navigation file
                    gt = np.array([float(x) for x in track_f.readline().split()])
                    try:
                        est = np.array([float(x) for x in sim_output_f.readline().split()])
                        est_cov = np.array([[float(x) for x in sim_output_f.readline().split()] for _ in range(dimension)])

                    # Empty file or not enough lines, the rest of the estimates are marked as Failed
                    except Exception as e:
                        print("Could not interpret estimate from output '%s'" % sim_output_f)
                        errors.append('Failed')
                        continue

                    # May happen the first line is read but no values, mark it is Failed as well
                    if len(est) == 0:
                        errors.append('Failed')
                        continue

                    # Compute and save the error
                    e = sum((gt - est)**2)
                    errors.append(e)

        # Save all sim errors
        all_errors.append(errors)

    # Compute the mean error at each timestep
    mean_errors = []
    for t in range(timesteps):

        # Handle failed values by ignoring them in the mean computation
        vals = [e[t] for e in all_errors if type(e[t]) != str]

        # If all Failed then mean Failed as well
        if len(vals) == 0:
            mean = 'Failed'
        else:
            mean = np.mean(vals)
        mean_errors.append(mean)

    assert(len(all_errors) == num_simulations)
    assert(len(mean_errors) == timesteps)

    # Save individual error files for each simulation
    for i in range(1, num_simulations+1):
        with open(error_output_filepath_base % i, 'w') as error_f:
            error_f.write('\n'.join(str(x) for x in all_errors[i-1]))

    # Save the mean error of all the simulations
    with open(mean_error_output_filepath, 'w') as mean_error_f:
        mean_error_f.write('\n'.join(str(x) for x in mean_errors))

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return



if __name__ == '__main__':
    create_sim_error_files(TRACK_FILEPATH_DEFAULT, SIM_OUTPUT_FILEPATH_BASE_DEFAULT, ERROR_OUTPUT_FILEPATH_BASE_DEFAULT, MEAN_ERROR_OUTPUT_FILEPATH_DEFAULT, NUM_SIMULATIONS_DEFAULT)