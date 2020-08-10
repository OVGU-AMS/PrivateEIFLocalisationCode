"""

"""

import os
import numpy as np

# Defaults for when this file is run
TRACK_FILEPATH_DEFAULT = "input/debug_track1.txt"
SIM_OUTPUT_FILEPATH_BASE_DEFAULT = "output/debug__nav_%03d.txt"
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

    all_errors = []
    # Get and save all the errors of each simulation to the track file
    for i in range(1, num_simulations+1):
        errors = []
        to_skip = False
        with open(track_filepath, 'r') as track_f:
            with open(sim_output_filepath_base % i, 'r') as sim_output_f:

                timesteps = int(track_f.readline())
                dimension = int(track_f.readline())
                init_state = np.array([float(x) for x in track_f.readline().split()])
                init_cov = np.array([[float(x) for x in track_f.readline().split()] for _ in range(dimension)])

                for t in range(timesteps):
                    gt = np.array([float(x) for x in track_f.readline().split()])
                    try:
                        est = np.array([float(x) for x in sim_output_f.readline().split()])
                        est_cov = np.array([[float(x) for x in sim_output_f.readline().split()] for _ in range(dimension)])
                    except Exception as e:
                        print("Could not interpret simulation output '%s'" % sim_output_f)
                        to_skip = True
                        break
                    
                    if np.allclose(est_cov, np.array([[0 for _ in range(dimension)] for _ in range(dimension)])):
                        print("Covariance too close to zero! Assuming fault and skipping simulation '%s'" % sim_output_f)
                        to_skip = True
                        break

                    if not all([e > 0 for e in np.linalg.eigvals(est_cov)]):
                        print("Covariance has non-positive eigenvalues! Assuming fault and skipping simulation '%s'" % sim_output_f)
                        to_skip = True
                        break

                    # Compute and save the error
                    e = np.linalg.norm(gt - est)
                    errors.append(e)

        if not to_skip:
            all_errors.append(errors)

    # Compute the mean error at each timestep
    mean_errors = []
    for t in range(timesteps):
        mean = np.mean([e[t] for e in all_errors])
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