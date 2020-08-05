"""

"""

import numpy as np


TRACK_FILEPATH = "../input/track1.txt"
SENSOR_FILEPATH_BASE = "../input/encoding_sim_%03d_sensor%d.txt"
NUMBER_OF_SIMS = 100

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

NUM_SENSORS = 4
FIRST_SENSOR_INDEX = 0

SENSOR_LOCATIONS = SENSOR_LOCATIONS[FIRST_SENSOR_INDEX:]


def generate_sim_inputs():
    for s in range(1, NUMBER_OF_SIMS+1):
        track_f = open(TRACK_FILEPATH, 'r')
        sensor_fs = [open(SENSOR_FILEPATH_BASE % (s, i), 'w') for i in range(1, NUM_SENSORS+1)]
        all_sensor_measurements = [[] for i in range(NUM_SENSORS)]

        timesteps = int(track_f.readline())
        dimension = int(track_f.readline())
        init_state = np.array([float(x) for x in track_f.readline().split()])
        init_cov = np.array([[float(x) for x in track_f.readline().split()] for _ in range(dimension)])

        for i in range(NUM_SENSORS):
            all_sensor_measurements[i].append("%d" % timesteps)
            all_sensor_measurements[i].append("%d" % dimension)
            all_sensor_measurements[i].append("%lf %lf" % (SENSOR_LOCATIONS[i][0], SENSOR_LOCATIONS[i][1]))

        for t in range(timesteps):
            gt = np.array([float(x) for x in track_f.readline().split()])
            for i in range(NUM_SENSORS):
                m = h(gt, SENSOR_LOCATIONS[i]) + np.random.normal(0, np.sqrt(SENSOR_VARIANCE))
                all_sensor_measurements[i].append("%lf" % m)


        track_f.close()
        for i in range(NUM_SENSORS):
            sensor_fs[i].write('\n'.join(all_sensor_measurements[i]))
            sensor_fs[i].close()
    return


def h(x, s):
    return np.linalg.norm(np.array([x[0], x[2]]) - s)


if __name__ == "__main__":
    generate_sim_inputs()