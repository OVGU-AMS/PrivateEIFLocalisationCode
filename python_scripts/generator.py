"""

"""

import os
import numpy as np


# These may remain constant, use num_sensors and first_sensor_index to select the sensor subset
SENSOR_VARIANCE = 5.0
SENSOR_LOCATIONS = [np.array([5.0, 5.0]), # Normal
                    np.array([40.0, 5.0]), 
                    np.array([5.0, 40.0]), 
                    np.array([40.0, 40.0]),
                    np.array([-30.0, -30.0]), # Big
                    np.array([75.0, -30.0]), 
                    np.array([-30.0, 75.0]), 
                    np.array([75.0, 75.0]),
                    np.array([-65.0, -65.0]), # Very big 
                    np.array([110.0, -65.0]), 
                    np.array([-65.0, 110.0]), 
                    np.array([110.0, 110.0]),
                    np.array([22.0, 22.0]), # Small
                    np.array([23.0, 22.0]), 
                    np.array([22.0, 23.0]), 
                    np.array([23.0, 23.0])]

# Defaults when this file run
TRACK_FILEPATH_DEFAULT = "input/debug_track1.txt"
SENSOR_FILEPATH_BASE_DEFAULT = "input/debug_sim_%03d_sensor%d.txt"
NUM_SIMS_DEFAULT = 1
FIRST_SENSOR_INDEX_DEFAULT = 0
NUM_SENSORS_DEFAULT = 4


# generate inputs for a particular simulation setup
def generate_sim_inputs(track_filepath, sensor_filepath_base, number_of_sims, first_sensor_index, number_of_sensors):
    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    for s in range(1, number_of_sims+1):
        track_f = open(track_filepath, 'r')
        sensor_fs = [open(sensor_filepath_base % (s, str(i)), 'w') for i in range(1, number_of_sensors+1)]
        all_sensor_measurements = [[] for i in range(number_of_sensors)]

        timesteps = int(track_f.readline())
        dimension = int(track_f.readline())
        init_state = np.array([float(x) for x in track_f.readline().split()])
        init_cov = np.array([[float(x) for x in track_f.readline().split()] for _ in range(dimension)])

        for i in range(number_of_sensors):
            all_sensor_measurements[i].append("%d" % timesteps)
            all_sensor_measurements[i].append("%d" % dimension)
            all_sensor_measurements[i].append("%lf %lf" % (SENSOR_LOCATIONS[first_sensor_index+i][0], SENSOR_LOCATIONS[first_sensor_index+i][1]))

        for t in range(timesteps):
            gt = np.array([float(x) for x in track_f.readline().split()])
            for i in range(number_of_sensors):
                m = h(gt, SENSOR_LOCATIONS[first_sensor_index+i]) + np.random.normal(0, np.sqrt(SENSOR_VARIANCE))
                # Ensure distance measurements are non-negative
                if m < 0:
                    m = 0
                all_sensor_measurements[i].append("%lf" % m)


        track_f.close()
        for i in range(number_of_sensors):
            sensor_fs[i].write('\n'.join(all_sensor_measurements[i]))
            sensor_fs[i].close()

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return


# Measurement function
def h(x, s):
    return np.linalg.norm(np.array([x[0], x[2]]) - s)


# On running this file, run with default params above
if __name__ == "__main__":
    generate_sim_inputs(TRACK_FILEPATH_DEFAULT, SENSOR_FILEPATH_BASE_DEFAULT, NUM_SIMS_DEFAULT, FIRST_SENSOR_INDEX_DEFAULT, NUM_SENSORS_DEFAULT)