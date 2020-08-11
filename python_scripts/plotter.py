import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Circle
import numpy as np
import os


# Single plot defaults when this file is run
TRACK_FILEPATH_DEFAULT = "../input/debug_track1.txt"
OUTPUT_FILEPATH_DEFAULT = "../output/debug_nav_001.txt"
SENSOR_LOCATIONS_DEFAULT = [np.array([5.0, 5.0]), # Normal
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
NUM_SENSORS_DEFAULT = 16

# 8888888          d8b 888
#   888            Y8P 888
#   888                888
#   888   88888b.  888 888888
#   888   888 "88b 888 888
#   888   888  888 888 888
#   888   888  888 888 Y88b.
# 8888888 888  888 888  "Y888




def init_matplotlib_params(save_not_show_fig, show_latex_fig):
    if save_not_show_fig:
        matplotlib.use("pgf")
        matplotlib.rcParams.update({
            "pgf.texsystem": "pdflatex",
            'font.family': 'serif',
            'text.usetex': True,
            'pgf.rcfonts': False,
        })
    elif show_latex_fig:
        matplotlib.rcParams.update({
            "pgf.texsystem": "pdflatex",
            'font.family': 'serif',
            'text.usetex': True,
            'pgf.rcfonts': False,
        })
    return



# 88888888888 d8b               d8b
#     888     Y8P               Y8P
#     888
#     888     888 88888b.d88b.  888 88888b.   .d88b.
#     888     888 888 "888 "88b 888 888 "88b d88P"88b
#     888     888 888  888  888 888 888  888 888  888
#     888     888 888  888  888 888 888  888 Y88b 888
#     888     888 888  888  888 888 888  888  "Y88888
#                                                 888
#                                            Y8b d88P
#                                             "Y88P"

def create_timing_plots(output_times_filepath_base, sensor_count_list, bitsize_list, width, height):

    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    # Reverse sort the bitzise counts so that large labels come first
    bitsize_list.sort(reverse=True)

    fig = plt.figure()
    fig.set_size_inches(w=width, h=height)
    ax = fig.add_subplot(111)

    sensor_counts = {}
    for s in sensor_count_list:
        sensor_counts[s] = []
        for b in bitsize_list:
            with open(output_times_filepath_base % (b, s)) as timing_f:
                time_vals = [float(x.strip()) for x in timing_f.read().split() if x.strip() != 'Failed']
                if len(time_vals) == 0:
                    mean = -1
                else:
                    # print("%d Run times for %d sensors using %d bit encryption:" % (len(time_vals), s, b), time_vals)
                    mean = np.mean(time_vals)
                sensor_counts[s].append(mean)

        assert(len(sensor_counts[s]) == len(bitsize_list))

        # Plot for x axis as Paillier bitsize
        #ax.plot([512, 1024, 2048], sensor_counts[s], marker='.')

    for b in range(len(bitsize_list)):
        ax.plot(sensor_count_list, [sensor_counts[i][b] for i in sensor_count_list], marker='x', label='%d bit encryption' % bitsize_list[b])

    ax.set_xlabel('Number of sensors')
    ax.set_ylabel('Runtime ($s$)')
    ax.set_xticks(sensor_count_list)
    ax.legend(loc=1)

    if matplotlib.get_backend() == 'pgf':
        plt.savefig('pictures/timing.pdf')
    else:
        plt.show()

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return

# 8888888888                                888 d8b
# 888                                       888 Y8P
# 888                                       888
# 8888888    88888b.   .d8888b .d88b.   .d88888 888 88888b.   .d88b.
# 888        888 "88b d88P"   d88""88b d88" 888 888 888 "88b d88P"88b
# 888        888  888 888     888  888 888  888 888 888  888 888  888
# 888        888  888 Y88b.   Y88..88P Y88b 888 888 888  888 Y88b 888
# 8888888888 888  888  "Y8888P "Y88P"   "Y88888 888 888  888  "Y88888
#                                                                 888
#                                                            Y8b d88P
#                                                             "Y88P"

def create_encoding_plots(output_filepath_base, encoding_types_list, sim_timesteps, width, height):

    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    fig = plt.figure()
    fig.set_size_inches(w=width, h=height)
    ax = fig.add_subplot(111)

    encoding_errors = {}
    for encoding in encoding_types_list:
        with open(output_filepath_base % encoding) as encoding_f:
                mean_errors = [float(x.strip()) if x.strip() != 'Failed' else -1 for x in encoding_f.read().split()]

        encoding_errors[encoding] = mean_errors
        ax.plot([0.5*x for x in list(range(sim_timesteps))], mean_errors, label='Modulus: %d, Fractional bits: %d' % encoding)

    ax.set_xlabel('Simulation time ($s$)')
    ax.set_ylabel('Average Simulation Error (RMSE)')
    #ax.set_xticks(sensor_count_list)
    ax.legend()

    if matplotlib.get_backend() == 'pgf':
        plt.savefig('pictures/encoding.pdf')
    else:
        plt.show()

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return

# 888                                          888
# 888                                          888
# 888                                          888
# 888       8888b.  888  888  .d88b.  888  888 888888 .d8888b
# 888          "88b 888  888 d88""88b 888  888 888    88K
# 888      .d888888 888  888 888  888 888  888 888    "Y8888b.
# 888      888  888 Y88b 888 Y88..88P Y88b 888 Y88b.       X88
# 88888888 "Y888888  "Y88888  "Y88P"   "Y88888  "Y888  88888P'
#                        888
#                   Y8b d88P
#                    "Y88P"

def create_distance_plots(output_filepath_base, distance_layout_list, distance_layout_labels, sim_timesteps, width, height):

    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    fig = plt.figure()
    fig.set_size_inches(w=width, h=height)
    ax = fig.add_subplot(111)

    layout_errors = {}
    for i,layout in enumerate(distance_layout_list):
        with open(output_filepath_base % layout) as distance_f:
                mean_errors = [float(x.strip()) if x.strip() != 'Failed' else -1 for x in distance_f.read().split()]

        layout_errors[layout] = mean_errors
        ax.plot([0.5*x for x in list(range(sim_timesteps))], mean_errors, label='%s' % distance_layout_labels[i])

    ax.set_xlabel('Simulation time ($s$)')
    ax.set_ylabel('Average Simulation Error (RMSE)')
    #ax.set_xticks(sensor_count_list)
    ax.legend()

    if matplotlib.get_backend() == 'pgf':
        plt.savefig('pictures/layouts.pdf')
    else:
        plt.show()

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return

#  .d8888b.  d8b                   888                        888          888
# d88P  Y88b Y8P                   888                        888          888
# Y88b.                            888                        888          888
#  "Y888b.   888 88888b.   .d88b.  888  .d88b.       88888b.  888  .d88b.  888888
#     "Y88b. 888 888 "88b d88P"88b 888 d8P  Y8b      888 "88b 888 d88""88b 888
#       "888 888 888  888 888  888 888 88888888      888  888 888 888  888 888
# Y88b  d88P 888 888  888 Y88b 888 888 Y8b.          888 d88P 888 Y88..88P Y88b.
#  "Y8888P"  888 888  888  "Y88888 888  "Y8888       88888P"  888  "Y88P"   "Y888
#                              888                   888
#                         Y8b d88P                   888
#                          "Y88P"                    888

def plot_sim(track_filepath, output_filepath, sensor_locations=None, num_sensors=None):
    fig = plt.figure()
    fig.set_size_inches(w=8, h=8)
    ax = fig.add_subplot(111)

    gts = []
    states = []
    covariances = []

    with open(track_filepath) as in_f:
        with open(output_filepath) as out_f:
            t = int(in_f.readline())
            dim = int(in_f.readline())
            init_state = np.array([float(x) for x in in_f.readline().split()])
            init_cov = np.array([[float(x) for x in in_f.readline().split()] for _ in range(dim)])

            for _ in range(t):
                gts.append(np.array([float(x) for x in in_f.readline().split()]))
                states.append(np.array([float(x) for x in out_f.readline().split()]))
                covariances.append(np.array([[float(x) for x in out_f.readline().split()] for _ in range(dim)]))


    # Plot initial estimate
    ax.scatter(init_state[0], init_state[2], marker='x', color='red')
    ax.add_artist(get_cov_ellipse(np.array([[init_cov[0][0], init_cov[0][2]],[init_cov[2][0], init_cov[2][2]]]), 
                                        np.array([init_state[0],init_state[2]]), 
                                        2, fill=False, linestyle='-', edgecolor='orange', zorder=1))


    # Plot ground truth
    ax.plot([x[0] for x in gts], [x[2] for x in gts], marker='.', color='gray')
    ax.plot([x[0] for x in states], [x[2] for x in states], marker='x', color='green')

    # plot filter estimates
    for state, cov in zip(states, covariances):
        ax.add_artist(get_cov_ellipse(np.array([[cov[0][0], cov[0][2]],[cov[2][0], cov[2][2]]]), 
                                        np.array([state[0],state[2]]), 
                                        2, fill=False, linestyle='-', edgecolor='royalblue', zorder=1))


    if sensor_locations != None:
        if num_sensors == None:
            num_sensors = len(sensor_locations)

        for s in range(num_sensors):
            ax.scatter(sensor_locations[s][0], sensor_locations[s][1], marker='o', color='red')

    # display the picture
    plt.show()


# 888    888          888
# 888    888          888
# 888    888          888
# 8888888888  .d88b.  888 88888b.   .d88b.  888d888 .d8888b
# 888    888 d8P  Y8b 888 888 "88b d8P  Y8b 888P"   88K
# 888    888 88888888 888 888  888 88888888 888     "Y8888b.
# 888    888 Y8b.     888 888 d88P Y8b.     888          X88
# 888    888  "Y8888  888 88888P"   "Y8888  888      88888P'
#                         888
#                         888
#                         888

# Plotting helpers, from https://scipython.com/book/chapter-7-matplotlib/examples/bmi-data-with-confidence-ellipses/
def get_cov_ellipse(cov, centre, nstd, **kwargs):
    """
    Return a matplotlib Ellipse patch representing the covariance matrix
    cov centred at centre and scaled by the factor nstd.

    """

    # Find and sort eigenvalues and eigenvectors into descending order
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # The anti-clockwise angle to rotate our ellipse by 
    vx, vy = eigvecs[:,0][0], eigvecs[:,0][1]
    theta = np.arctan2(vy, vx)

    # Width and height of ellipse to draw
    width, height = 2 * nstd * np.sqrt(eigvals) # eigvals positive because covariance is positive semi definite
    return Ellipse(xy=centre, width=width, height=height,
                   angle=np.degrees(theta), **kwargs)



if __name__ == '__main__':
    init_matplotlib_params(False, False)
    plot_sim(TRACK_FILEPATH_DEFAULT, OUTPUT_FILEPATH_DEFAULT, SENSOR_LOCATIONS_DEFAULT, NUM_SENSORS_DEFAULT)

    #create_timing_plots('output/timing_%d_%d_nav_times.txt', [2,3,4,5], [512, 1024, 1536, 2048, 2560])
    #create_encoding_plots('output_evaluation/encoding_%d_%d_nav_mean_errors.txt', [(64, 16), (128, 32), (256, 64)], 50)
    #create_distance_plots('output_evaluation/layout_%s_nav_mean_errors.txt', ['small', 'normal', 'big', 'verybig'], ['Small', 'Normal', 'Large', 'Very Large'], 50)
