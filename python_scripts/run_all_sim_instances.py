"""

"""

import subprocess
import numpy as np


NUM_SENSORS = 4
NUM_SIMULATIONS = 1

INPUT_TRACK_FILEPATH = "input/track1.txt"

OUTPUT_FILEPATH_BASE = "output/encoding_nav%03d.txt"
# OUTPUT_FILEPATH_BASE = "../output/timing_nav%03d_%dsensors.txt"
# OUTPUT_FILEPATH_BASE = "../output/layout_small_nav%03d.txt"
# OUTPUT_FILEPATH_BASE = "../output/layout_normal_nav%03d.txt"
# OUTPUT_FILEPATH_BASE = "../output/layout_big_nav%03d.txt"
# OUTPUT_FILEPATH_BASE = "../output/layout_verybig_nav%03d.txt"

SENSOR_FILEPATH_BASE = "input/encoding_sim_%03d_sensor%s.txt"
#SENSOR_FILEPATH_BASE = "../input/timing_sim_%03d_sensor%d.txt"
#SENSOR_FILEPATH_BASE = "../input/layout_small_sim_%03d_sensor%d.txt"
#SENSOR_FILEPATH_BASE = "../input/layout_normal_sim_%03d_sensor%d.txt"
#SENSOR_FILEPATH_BASE = "../input/layout_big_sim_%03d_sensor%d.txt"
#SENSOR_FILEPATH_BASE = "../input/layout_verybig_sim_%03d_sensor%d.txt"


def run_all_simulation_repeats(num_simulations, num_sensors, input_track_fp, output_fp_base, sensor_fp_base):
	cwd = '../'
	for i in range(1, num_simulations+1):
		print("Running simulation %d" % i)
		args = ['mpirun', '-np', str(num_sensors+1), 'build/sim', input_track_fp, output_fp_base % i, sensor_fp_base % (i, "%d")]
		print(args)

		# shell=False (do not use the the shell -> no piping, redirecting, etc, but much faster)
		# text=True (everything returned and output by the subprocess are strings)

		p = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False, cwd=cwd)

		print(p.stdout)
		print(p.stderr)


	return




if __name__ == '__main__':
	run_all_simulation_repeats(NUM_SIMULATIONS, NUM_SENSORS, INPUT_TRACK_FILEPATH, OUTPUT_FILEPATH_BASE, SENSOR_FILEPATH_BASE)