import numpy as np
import os


TRACK_FILEPATH_DEFAULT = "input/track_%03d.txt"
NUM_SIMS_DEFAULT = 50


# Create input tracks to use with simulation
def generate_input_tracks(track_filepath, number_of_sims):
    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    # System variables
    q = 0.01 # Noise strength
    t = 0.5 # Time step
    F = np.array([[1, t, 0, 0],[0, 1, 0, 0],[0, 0, 1, t],[0, 0, 0, 1]])
    Q = q*np.array([[t**3/3,t**2/2,0,0],[t**2/2,t,0,0],[0,0,t**3/3,t**2/2],[0,0,t**2/2,t]])

    # Generation variables
    timesteps = 50
    dimension = 4
    init_state = np.array([0, 1, 0, 1])
    init_cov = np.array([[100, 0, 0, 0],
                         [0, 100, 0, 0],
                         [0, 0, 100, 0],
                         [0, 0, 0, 100]])

    for sim in range(1, number_of_sims+1):
        with open(track_filepath % sim, 'w') as output_f:

            output_f.write("%d\n" % timesteps)
            output_f.write("%d\n" % dimension)
            output_f.write("%lf %lf %lf %lf\n" % (init_state[0], init_state[1], init_state[2], init_state[3]))
            for r in range(dimension):
                    output_f.write("%lf %lf %lf %lf\n" % (init_cov[r][0], init_cov[r][1], init_cov[r][2], init_cov[r][3]))

            state = init_state.copy()
            cov = init_cov.copy()
            for t in range(timesteps):
                # Noise
                w = np.random.multivariate_normal(np.array([0,0,0,0]), Q)
                # Next state
                state = F @ state + w
                # write to file
                output_f.write("%lf %lf %lf %lf\n" % (state[0], state[1], state[2], state[3]))

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return




if __name__ == '__main__':
    generate_input_tracks(TRACK_FILEPATH_DEFAULT, NUM_SIMS_DEFAULT)