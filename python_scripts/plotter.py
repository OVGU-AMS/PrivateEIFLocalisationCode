import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Circle
import numpy as np
import os


# Single plot defaults when this file is run
TRACK_FILEPATH_DEFAULT = "input/track_117.txt"
OUTPUT_FILEPATH_DEFAULT = "output/layout_big_nav_117.txt"
SENSOR_FILEPATH_BASE = "input/layout_big_sim_117_sensor%d.txt"
NUM_SENSORS_DEFAULT = 4
ADDITIONAL_SENSOR_LOCATIONS_DEFAULT = []


# 8888888          d8b 888
#   888            Y8P 888
#   888                888
#   888   88888b.  888 888888
#   888   888 "88b 888 888
#   888   888  888 888 888
#   888   888  888 888 Y88b.
# 8888888 888  888 888  "Y888




def init_matplotlib_params(save_not_show_fig, show_latex_fig):
    fontsize = 9
    linewidth = 1.0
    gridlinewidth = 0.7

    if save_not_show_fig:
        matplotlib.use("pgf")

    if save_not_show_fig or show_latex_fig:
        matplotlib.rcParams.update({
            "pgf.texsystem": "pdflatex",   
            'font.family': 'serif',         # Use serif/main font for text elements
            'text.usetex': True,            # Use inline maths for ticks
            'pgf.rcfonts': False,           # Don't setup fonts from matplotlib rc params
            'font.size': fontsize,
            'axes.titlesize': fontsize,
            'axes.labelsize': fontsize,
            'xtick.labelsize': fontsize,
            'ytick.labelsize': fontsize,
            'legend.fontsize': fontsize,
            'figure.titlesize': fontsize,
            # line width
            'lines.linewidth': linewidth,
            'grid.linewidth': gridlinewidth
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

def create_timing_plots(output_times_filepath_base, sensor_count_list, bitsize_list, width):
    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    # Reverse sort the bitzise counts so that large labels come first
    bitsize_list.sort(reverse=True)

    # Fixed width for template, the rest if eye-balled
    fig = plt.figure()
    fig.set_size_inches(w=width, h=2)
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.175, right=0.825, bottom=0.16, top=0.63)
    ax.grid(linestyle='dashed')
    ax.set_axisbelow(True)

    sensor_counts = {}
    for s in sensor_count_list:
        sensor_counts[s] = []
        for b in bitsize_list:
            with open(output_times_filepath_base % (b, s)) as timing_f:
                time_vals = [float(x.strip()) for x in timing_f.read().split() if x.strip() != 'Failed']
                if len(time_vals) == 0:
                    mean = -1
                else:
                    mean = np.mean(time_vals)
                sensor_counts[s].append(mean)

        assert(len(sensor_counts[s]) == len(bitsize_list))

    plot_handles = []
    for b in range(len(bitsize_list)):
        ph, = ax.plot(sensor_count_list, [sensor_counts[i][b] for i in sensor_count_list], marker='x', label=r'%d' % bitsize_list[b])
        plot_handles.append(ph)

    # Set the ticks to be all of the sensor amounts tested
    ax.set_xticks(sensor_count_list)

    # Used gloabl axis labels to match the other figures
    fig.supxlabel(r'Number of sensors')
    fig.supylabel(r'Runtime ($s$)')

    # Legend
    fig.legend(handles=plot_handles, title='Key Length (bits)', loc='upper center', ncol=2)

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

def create_encoding_plots(output_filepath_base, encoding_types_list, sim_timesteps, width, eif_filepath):

    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    fig = plt.figure()
    fig.set_size_inches(w=2*width, h=2*width)
    ax = fig.add_subplot(111)

    # Plot the private EIF with varying encodings
    encoding_errors = {}
    for encoding in encoding_types_list:
        with open(output_filepath_base % encoding) as encoding_f:
                mean_errors = [float(x.strip()) if x.strip() != 'Failed' else -1 for x in encoding_f.read().split()]

        encoding_errors[encoding] = mean_errors
        ax.plot([0.5*x for x in list(range(sim_timesteps))], mean_errors, label=r'Fractional bits: %d' % encoding)

    # Plot the normal EIF
    if eif_filepath != None:
        with open(eif_filepath) as eif_encoding_f:
                mean_errors = [float(x.strip()) if x.strip() != 'Failed' else -1 for x in eif_encoding_f.read().split()]
        ax.plot([0.5*x for x in list(range(sim_timesteps))], mean_errors, label=r'EIF')

    ax.set_xlabel(r'Simulation time ($s$)')
    ax.set_ylabel(r'Average Simulation Error ($m$)')
    ax.legend()

    if matplotlib.get_backend() == 'pgf':
        plt.savefig('pictures/encoding_errors.pdf')
    else:
        plt.show()

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return

# 8888888b.  d8b          888
# 888  "Y88b Y8P          888
# 888    888              888
# 888    888 888 .d8888b  888888  8888b.  88888b.   .d8888b .d88b.  .d8888b
# 888    888 888 88K      888        "88b 888 "88b d88P"   d8P  Y8b 88K
# 888    888 888 "Y8888b. 888    .d888888 888  888 888     88888888 "Y8888b.
# 888  .d88P 888      X88 Y88b.  888  888 888  888 Y88b.   Y8b.          X88
# 8888888P"  888  88888P'  "Y888 "Y888888 888  888  "Y8888P "Y8888   88888P'




def create_distance_plots(output_filepath_base, distance_layout_list, distance_layout_labels, sim_timesteps, width, eif_filepath_base):
    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    # Width fixed for template, the rest is eye-balled
    fig, axs = plt.subplots(4,1, figsize=(width, 4.5), sharex=True, sharey=False)
    plt.subplots_adjust(wspace=0, hspace=0.375, top=0.81, bottom=0.0825, left=0.16, right=0.84)
    
    enc_plot_handles = []
    eif_plot_handles = []
    diff_plot_handles = []

    # Axis for each layout
    for i,layout in enumerate(distance_layout_list):
        ax = axs.flat[i]
        ax.set_title(distance_layout_labels[i])
        ax.grid(linestyle='dashed')
        if i < len(distance_layout_list)-1:
            ax.tick_params(labelcolor='none', bottom=False)

        # Plot our filter errors
        with open(output_filepath_base % layout) as distance_f:
            mean_enc_errors = [float(x.strip()) if x.strip() != 'Failed' else -1 for x in distance_f.read().split()]
        # Note [1:] skips the first estimate as it is the initial state and same for all
        ph_enc, = ax.plot([x for x in list(range(sim_timesteps))][1:], mean_enc_errors[1:])
        enc_plot_handles.append(ph_enc)

        # If a standard EIF errors file is give, plot that and the difference between it and ours as well
        if eif_filepath_base != None:
            with open(eif_filepath_base % layout) as eif_distance_f:
                mean_eif_errors = [float(x.strip()) if x.strip() != 'Failed' else -1 for x in eif_distance_f.read().split()]
            mean_diffs = [abs(x-y) for x,y in zip(mean_enc_errors, mean_eif_errors)]
            # Note [1:] skips the first estimate as it is the initial state and same for all
            ph_eif, = ax.plot([x for x in list(range(sim_timesteps))][1:], mean_eif_errors[1:], linestyle='--')
            ph_dif, = ax.plot([x for x in list(range(sim_timesteps))][1:], mean_diffs[1:], linestyle='--')
            eif_plot_handles.append(ph_eif)
            diff_plot_handles.append(ph_dif)

    # Shared axis labels
    fig.supxlabel(r'Filter Timesteps')
    fig.supylabel(r'Simulation Error (MSE)')

    # Legend only uses lines from first plot (all are the same colours)
    fig.legend((enc_plot_handles[0], eif_plot_handles[0], diff_plot_handles[0]), (r'Privacy-preserving Filter', r'Standard EIF', r'Difference'), loc='upper center', ncol=1)

    if matplotlib.get_backend() == 'pgf':
        plt.savefig('pictures/layout_errors.pdf')
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

def plot_layouts_and_track(track_filepath, sensor_lists, sensor_list_labels):

    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

    # Get track file measurements
    with open(track_filepath % 1) as track_f:
        ground_truth = []
        timesteps = int(track_f.readline())
        dimenstions = int(track_f.readline())
        init_state = np.array([float(x) for x in track_f.readline().split()])
        init_cov = np.array([[float(x) for x in track_f.readline().split()] for _ in range(dimenstions)])
        for _ in range(timesteps):
            ground_truth.append(np.array([float(x) for x in track_f.readline().split()]))

    # Fixed width from template
    width = 2.95

    # Make subplots. Special sizes to accomodate Automatica column width
    fig, axs = plt.subplots(2,2, figsize=(width, 3), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=-0.25, hspace=0.35, top=0.71, bottom=0.13, left=0.1, right=0.9)

    plots = []
    scatters = []
    for i,ax in enumerate(axs.flat):

        ax.set_title(sensor_list_labels[i])
        ax.set_aspect(aspect='equal')

        # Plot sensor positions and groundtruth
        s = ax.scatter([x[0] for x in sensor_lists[i]], [x[1] for x in sensor_lists[i]], marker='.', color='red')
        p, = ax.plot([x[0] for x in ground_truth], [x[2] for x in ground_truth])
        i_s = ax.scatter([init_state[0]], [init_state[2]], marker='.', color='orange')
        ax.add_artist(get_cov_ellipse(np.array([[init_cov[0][0], init_cov[0][2]],[init_cov[2][0], init_cov[2][2]]]), 
                                        np.array([init_state[0],init_state[2]]), 
                                        2, fill=False, linestyle='-', edgecolor='orange', zorder=1))

        # Save the plots for a shared legend
        scatters.append(s)
        plots.append(p)

    # Legend
    fig.legend((plots[0], i_s, scatters[0]), (r'Ground Truth', r'Initial Estimate and Covariance', r'Sensors'), loc='upper center', ncol=1)

    # Shared axis labels
    fig.supxlabel(r'Location $x$')   
    fig.supylabel(r'Location $y$')

    # Hide ticks from intermediate axes
    for a in [axs[0][0], axs[0][1]]:
        a.tick_params(labelcolor='none', bottom=False)

    # Hide ticks from intermediate axes
    for a in [axs[0][1], axs[1][1]]:
        a.tick_params(labelcolor='none', left=False)

    # Save
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

def plot_sim(track_filepath, output_filepath, sensor_filepath_base=None, num_sensors=None, additional_sensor_locations=None):
    
    # Always run from top project folder, if currently in python folder, move up
    dir_moved = False
    if os.getcwd().endswith('python_scripts'):
        os.chdir('../')
        dir_moved = True

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

    plot_sensors = False
    if sensor_filepath_base != None:
        plot_sensors = True
        locs = []
        for i in range(1, num_sensors+1):
            with open(sensor_filepath_base % i) as sen_f:
                t = int(sen_f.readline())
                sen_f.readline()
                loc = np.array([float(x.strip()) for x in sen_f.readline().split()])
                locs.append(loc)

    if additional_sensor_locations != None:
        plot_sensors = True
        for extra_loc in additional_sensor_locations:
            locs.append(extra_loc)

    if plot_sensors:
        ax.scatter([l[0] for l in locs], [l[1] for l in locs], marker='o', color='red')

    # display the picture
    plt.show()

    # If directory was moved up, move it back down before ending
    if dir_moved:
        os.chdir('./python_scripts')

    return


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
    plot_sim(TRACK_FILEPATH_DEFAULT, OUTPUT_FILEPATH_DEFAULT, SENSOR_FILEPATH_BASE, NUM_SENSORS_DEFAULT, ADDITIONAL_SENSOR_LOCATIONS_DEFAULT)

    #create_timing_plots('output/timing_%d_%d_nav_times.txt', [2,3,4,5], [512, 1024, 1536, 2048, 2560])
    #create_encoding_plots('output_evaluation/encoding_%d_%d_nav_mean_errors.txt', [(64, 16), (128, 32), (256, 64)], 50)
    #create_distance_plots('output_evaluation/layout_%s_nav_mean_errors.txt', ['small', 'normal', 'big', 'verybig'], ['Small', 'Normal', 'Large', 'Very Large'], 50)
    # sensor_locations = [[np.array([5.0, 5.0]), # Normal
    #                     np.array([40.0, 5.0]), 
    #                     np.array([5.0, 40.0]), 
    #                     np.array([40.0, 40.0])],
    #                     [np.array([-30.0, -30.0]), # Big
    #                     np.array([75.0, -30.0]), 
    #                     np.array([-30.0, 75.0]), 
    #                     np.array([75.0, 75.0])],
    #                     [np.array([-65.0, -65.0]), # Very big 
    #                     np.array([110.0, -65.0]), 
    #                     np.array([-65.0, 110.0]), 
    #                     np.array([110.0, 110.0])],
    #                     [np.array([-100.0, -100.0]), # Huge
    #                     np.array([145.0, -100.0]), 
    #                     np.array([-100.0, 145.0]), 
    #                     np.array([145.0, 145.0])], 
    #                     [np.array([22.0, 22.0]), # Small                        
    #                     np.array([23.0, 22.0]), 
    #                     np.array([22.0, 23.0]), 
    #                     np.array([23.0, 23.0])]]
    # plot_layouts_and_track('input/debug_track_1.txt', sensor_locations, ['Normal', 'Big', 'Quite Big', 'Very Big', 'Small'])
