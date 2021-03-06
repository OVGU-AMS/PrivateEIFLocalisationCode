import numpy as np
import os


TRACK_FILEPATH_DEFAULT = "input/debug_track_%03d.txt"
SENSOR_FILEPATH_BASE_DEFAULT = "input/encoding_128_32_sim_%03d_sensor%d.txt"
OUTPUT_FILEPATH_BASE = "output/eif_encoding_128_32_nav_%03d.txt"
NUM_SIMS_DEFAULT = 50
NUM_SENSORS_DEFAULT = 4


# Run a non encryption extended information filter
def run_normal_eif(track_filepath, sensor_filepath_base, output_filepath_base, number_of_sims, number_of_sensors, sim_range_to_run=None):
    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    # Filter variables
    q = 0.01 # Noise strength
    t = 0.5 # Time step
    F = np.array([[1, t, 0, 0],[0, 1, 0, 0],[0, 0, 1, t],[0, 0, 0, 1]])
    Q = q*np.array([[t**3/3,t**2/2,0,0],[t**2/2,t,0,0],[0,0,t**3/3,t**2/2],[0,0,t**2/2,t]])
    R = 5

    # For manual running control, allow choosing exactly which sim numbers to run (in case of pc crashed/restarts/updates)
    start_sim = 1
    end_sim = number_of_sims+1
    if sim_range_to_run:
        start_sim = sim_range_to_run[0]
        end_sim = sim_range_to_run[1]

    for sim in range(start_sim, end_sim):
        # Open track file
        track_f = open(track_filepath % sim, 'r')
        timesteps = int(track_f.readline())
        dimension = int(track_f.readline())
        init_state = np.array([float(x) for x in track_f.readline().split()])
        init_cov = np.array([[float(x) for x in track_f.readline().split()] for _ in range(dimension)])

        # Open all the sensor files
        sensor_files = [open(sensor_filepath_base % (sim, i)) for i in range(1, number_of_sensors+1)]
        sensor_locs = []

        for sen_f in sensor_files:
            t = int(sen_f.readline())
            d = int(sen_f.readline())
            sensor_locs.append(np.array([float(x.strip()) for x in sen_f.readline().split()]))

        with open(output_filepath_base % sim, 'w') as output_f:
            state = init_state.copy()
            cov = init_cov.copy()
            for t in range(timesteps):

                # Prediction
                state = F @ state
                cov = F @ cov @ F.T + Q

                # Sensors compute EIF vars
                hrhs = []
                hrzs = []
                for s in range(number_of_sensors):
                    z = float(sensor_files[s].readline())
                    nav_pos = np.array([state[0], state[2]])
                    sen_pos = sensor_locs[s]
                    H = np.array([(nav_pos[0] - sen_pos[0]) / np.sqrt((nav_pos[0] - sen_pos[0])**2 + (nav_pos[1] - sen_pos[1])**2), 
                                  0,
                                  (nav_pos[1] - sen_pos[1]) / np.sqrt((nav_pos[0] - sen_pos[0])**2 + (nav_pos[1] - sen_pos[1])**2),
                                  0])

                    hrhs.append(np.outer((H.T * (1/R)), H))
                    hrzs.append(H.T * (1/R) * (z - np.linalg.norm(nav_pos - sen_pos) + H @ state))

                # Update
                info_matrix = np.linalg.inv(cov)
                info_vector = info_matrix @ state

                info_matrix = info_matrix + sum(hrhs)
                info_vector = info_vector + sum(hrzs)

                cov = np.linalg.inv(info_matrix)
                state = cov @ info_vector

                output_f.write("%lf %lf %lf %lf\n" % (state[0], state[1], state[2], state[3]))
                for r in range(dimension):
                    output_f.write("%lf %lf %lf %lf\n" % (cov[r][0], cov[r][1], cov[r][2], cov[r][3]))
                
        # Close track and sensor files before next simulation
        track_f.close()
        for i in sensor_files:
            i.close()

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return




if __name__ == '__main__':
    run_normal_eif(TRACK_FILEPATH_DEFAULT, SENSOR_FILEPATH_BASE_DEFAULT, OUTPUT_FILEPATH_BASE, NUM_SIMS_DEFAULT, NUM_SENSORS_DEFAULT)