"""

"""

import os
import subprocess
import time
import numpy as np


# Defaults for when this file is run
TRACK_FILEPATH_DEFAULT = "input/debug_track_%03d.txt"
SENSOR_FILEPATH_BASE_DEFAULT = "input/debug_sim_%03d_sensor%s.txt"
OUTPUT_FILEPATH_BASE_DEFAULT = "output/debug_nav_%03d.txt"
RUNTIME_OUTPUT_FILEPATH_DEFAULT = "output/debug_nav_times.txt"
NUM_SENSORS_DEFAULT = 4
PAILLIER_BITSIZE_DEFAULT = 1024
ENCODING_FRAC_BITSIZE_DEFAULT = 32
NUM_SIMULATIONS_DEFAULT = 1


def run_simulation_repeats(track_filepath, sensor_filepath_base, output_filepath_base, runtimes_filepath, num_sensors, paillier_bitsize, encoding_frac_bitsize, num_simulations):
    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    runtimes = []
    for i in range(1, num_simulations+1):
        # Filepath commanline args for the simulation
        track_fp = track_filepath % i
        out_fp = output_filepath_base % i
        sensors_fpb = sensor_filepath_base % (i, "%d")

        # Command
        args = ['mpirun', '-np', str(num_sensors+1), 'build/sim', str(paillier_bitsize), str(encoding_frac_bitsize), out_fp, track_fp, sensors_fpb]
        print("Running simulation %d:" % i, ' '.join(args))

        # Run simulation
        p = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)

        # Get runtime from program stdout
        try:
            runtimes.append(float(p.stdout))
        except:
            print("Failed to convert whole output to float!")
            try:
                runtimes.append(float(p.stdout.strip().split('\n')[-1]))
            except:
                print("Failed to convert end of output to float as well!")
                print('stderr:', p.stderr)
                runtimes.append('Failed')

        print('stdout:', p.stdout)

    # Write all simulation runtimes to time_output file
    with open(runtimes_filepath, 'w') as runtimes_f:
        runtimes_f.write('\n'.join([str(t) for t in runtimes]))

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return




if __name__ == '__main__':
    run_simulation_repeats(TRACK_FILEPATH_DEFAULT, 
                            SENSOR_FILEPATH_BASE_DEFAULT, 
                            OUTPUT_FILEPATH_BASE_DEFAULT, 
                            RUNTIME_OUTPUT_FILEPATH_DEFAULT, 
                            NUM_SENSORS_DEFAULT,
                            PAILLIER_BITSIZE_DEFAULT, 
                            ENCODING_FRAC_BITSIZE_DEFAULT, 
                            NUM_SIMULATIONS_DEFAULT)