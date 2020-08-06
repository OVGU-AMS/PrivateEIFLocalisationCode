"""

"""

import subprocess
import numpy as np


# Defaults for when this file is run
TRACK_FILEPATH_DEFAULT = "input/debug_track1.txt"
SENSOR_FILEPATH_BASE_DEFAULT = "input/debug_sim_%03d_sensor%s.txt"
OUTPUT_FILEPATH_BASE_DEFAULT = "output/debug_nav%03d.txt"
RUNTIME_OUTPUT_FILEPATH_DEFAULT = "output/debug_nav_times.txt"
NUM_SENSORS_DEFAULT = 4
PAILLIER_BITSIZE_DEFAULT = 1024
ENCODING_MOD_BITSIZE_DEFAULT = 128
ENCODING_FRAC_BITSIZE_DEFAULT = 32
NUM_SIMULATIONS_DEFAULT = 1


def run_simulation_repeats(track_filepath, sensor_filepath_base, output_filepath_base, runtimes_filepath, num_sensors, paillier_bitsize, encoding_mod_bitsize, encoding_frac_bitsize, num_simulations):
	cwd = '../'
	runtimes = []
	for i in range(1, num_simulations+1):
		# Filepath commanline args for the simulation
		track_fp = track_filepath
		out_fp = output_filepath_base % i
		sensors_fpb = sensor_filepath_base % (i, "%d")

		# Command
		args = ['mpirun', '-np', str(num_sensors+1), 'build/sim', str(paillier_bitsize), str(encoding_mod_bitsize), str(encoding_frac_bitsize), out_fp, track_fp, sensors_fpb]
		print("Running simulation %d:" % i, ' '.join(args))

		# Run simulation
		# shell=False (do not use the the shell -> no piping, redirecting, etc, but much faster)
		p = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False, cwd=cwd)

		# Get runtime from program stdout
		try:
			runtimes.append(float(p.stdout))
		except:
			print("Failed to convert output to float!")

		print('stdout:', p.stdout)

	# Write all simulation runtimes to time_output file
	with open(cwd+runtimes_filepath, 'w') as runtimes_f:
		runtimes_f.write('\n'.join([str(t) for t in runtimes]))


	return




if __name__ == '__main__':
	run_simulation_repeats(TRACK_FILEPATH_DEFAULT, 
							SENSOR_FILEPATH_BASE_DEFAULT, 
							OUTPUT_FILEPATH_BASE_DEFAULT, 
							RUNTIME_OUTPUT_FILEPATH_DEFAULT, 
							NUM_SENSORS_DEFAULT,
							PAILLIER_BITSIZE_DEFAULT, 
							ENCODING_MOD_BITSIZE_DEFAULT, 
							ENCODING_FRAC_BITSIZE_DEFAULT, 
							NUM_SIMULATIONS_DEFAULT)