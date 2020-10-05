"""

"""

import numpy as np
import generator
import runner
import base_runner
import evaluator
import plotter

# Process (gen., sim., eval., etc.) repeats 
SIM_REPEATS = 5

# Which part of process to do
DO_MEASUREMENT_GEN = False
DO_SIM_RUN = False
DO_SIM_EVALUATE = False
DO_PLOT_CREATE = True

# Which scenarios to run
DO_ENCODING = False
DO_TIMING = True
DO_DISTANCE = True

# EIF base simulations
ENCODING_ONLY_EIF_BASE = False
DISTANCE_ONLY_EIF_BASE = False

# Plot generation and defaults
SAVE_NOT_SHOW_FIG = True
SHOW_LATEX_FIG = True
FIG_WIDTH_DEFAULT = 2.95

# Generate measurements, run simulations, evaluate results, and plot results (or any subset of these)
def run_all(sim_repeats, do_measurement_gen, do_sim_run, do_sim_evaluate, do_plot_create, do_encoding, do_timing, do_distance):

    """
    8888888b.
    888   Y88b
    888    888
    888   d88P 8888b.  888d888 8888b.  88888b.d88b.  .d8888b
    8888888P"     "88b 888P"      "88b 888 "888 "88b 88K
    888       .d888888 888    .d888888 888  888  888 "Y8888b.
    888       888  888 888    888  888 888  888  888      X88 d8b
    888       "Y888888 888    "Y888888 888  888  888  88888P' Y8P



    """

    # Shared params
    sensor_locations = [np.array([5.0, 5.0]), # Normal
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
                        np.array([-100.0, -100.0]), # Huge
                        np.array([145.0, -100.0]), 
                        np.array([-100.0, 145.0]), 
                        np.array([145.0, 145.0]), 
                        np.array([22.0, 22.0]), # Small
                        np.array([23.0, 22.0]), 
                        np.array([22.0, 23.0]), 
                        np.array([23.0, 23.0]),]

    # Encoding test params
    encodings_to_test = [8, 16, 32]
    encoding_plot_eif_base = True
    encoding_fig_width = FIG_WIDTH_DEFAULT

    # Timing test params
    timing_paillier_bitsizes_to_test = [512, 1024, 1536, 2048, 2560]
    min_timing_sensors_to_test = 2
    max_timing_sensors_to_test = 5
    timing_fig_width = FIG_WIDTH_DEFAULT

    # Distance test params
    layouts = ['normal', 'big', 'verybig', 'huge']
    layout_labels = ['Normal', 'Big', 'Quite Big', 'Very Big']
    layouts_plot_eif_base = True
    layouts_fig_width = FIG_WIDTH_DEFAULT

    """
    888b     d888                                                                                     888
    8888b   d8888                                                                                     888
    88888b.d88888                                                                                     888
    888Y88888P888  .d88b.   8888b.  .d8888b  888  888 888d888 .d88b.  88888b.d88b.   .d88b.  88888b.  888888       .d88b.   .d88b.  88888b.
    888 Y888P 888 d8P  Y8b     "88b 88K      888  888 888P"  d8P  Y8b 888 "888 "88b d8P  Y8b 888 "88b 888         d88P"88b d8P  Y8b 888 "88b
    888  Y8P  888 88888888 .d888888 "Y8888b. 888  888 888    88888888 888  888  888 88888888 888  888 888         888  888 88888888 888  888
    888   "   888 Y8b.     888  888      X88 Y88b 888 888    Y8b.     888  888  888 Y8b.     888  888 Y88b.       Y88b 888 Y8b.     888  888 d8b
    888       888  "Y8888  "Y888888  88888P'  "Y88888 888     "Y8888  888  888  888  "Y8888  888  888  "Y888       "Y88888  "Y8888  888  888 Y8P
                                                                                                                       888
                                                                                                                  Y8b d88P
                                                                                                                   "Y88P"
    """

    if do_measurement_gen:
        # Genereate measurements for encoding plots
        if do_encoding:
            print("Genereating measurements for encoding plot...")
            for frac_bits in encodings_to_test:
                generator.generate_sim_inputs("input/track1.txt", "input/encoding_" + str(frac_bits) + "_sim_%03d_sensor%s.txt", sim_repeats, sensor_locations, 0, 4)

        # Genereate measurements for timing plots
        if do_timing:
            print("Genereating measurements for timing plot...")
            generator.generate_sim_inputs("input/track1.txt", "input/timing_sim_%03d_sensor%s.txt", sim_repeats, sensor_locations, 0, 8)

        # Generate measurements for distance plots
        if do_distance:
            print("Genereating measurements for distance/layout plot...")
            for i,layout in enumerate(layouts):
                sensor_fpb = "input/layout_" + layout + "_sim_%03d_sensor%s.txt"
                sensor_start_index = i*4
                generator.generate_sim_inputs("input/track1.txt", sensor_fpb, sim_repeats, sensor_locations, sensor_start_index, 4)

        print("Finished genereating measurements.\n")


    """
    8888888b.                                  d8b
    888   Y88b                                 Y8P
    888    888
    888   d88P 888  888 88888b.       .d8888b  888 88888b.d88b.  .d8888b
    8888888P"  888  888 888 "88b      88K      888 888 "888 "88b 88K
    888 T88b   888  888 888  888      "Y8888b. 888 888  888  888 "Y8888b.
    888  T88b  Y88b 888 888  888           X88 888 888  888  888      X88 d8b
    888   T88b  "Y88888 888  888       88888P' 888 888  888  888  88888P' Y8P



    """

    if do_sim_run:
        # Run all encoding sims
        if do_encoding:
            print("Running simulations for encoding plot...")
            for frac_bits in encodings_to_test:

                in_fpb = "input/encoding_" + str(frac_bits) + "_sim_%03d_sensor%s.txt"
                out_fpb = "output/encoding_" + str(frac_bits) + "_nav_%03d.txt"
                out_times_fp = "output/encoding_" + str(frac_bits) + "_nav_times.txt"

                if not ENCODING_ONLY_EIF_BASE:
                    runner.run_simulation_repeats("input/track1.txt",
                                                    in_fpb,
                                                    out_fpb,
                                                    out_times_fp,
                                                    4, 1024, frac_bits, sim_repeats)

            # Run the normal EIF on only on eof the inputs, since encoding isn't used in the normal filter
            frac_bits = encodings_to_test[0]
            in_eif_fpb = "input/encoding_" + str(frac_bits) + "_sim_%03d_sensor%s.txt"
            out_eif_fpb = "output/eif_encoding_" + str(frac_bits) + "_nav_%03d.txt"
            base_runner.run_normal_eif("input/track1.txt", in_eif_fpb, out_eif_fpb, sim_repeats, 4)

        # Run all timing sims
        if do_timing:
            print("Running simulations for timing plot...")
            for bitsize in timing_paillier_bitsizes_to_test:
                for sensors in range(min_timing_sensors_to_test, max_timing_sensors_to_test+1):

                    out_fpb = "output/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_%03d.txt"
                    out_times_fp = "output/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_times.txt"
                    runner.run_simulation_repeats("input/track1.txt",
                                                    "input/timing_sim_%03d_sensor%s.txt",
                                                    out_fpb,
                                                    out_times_fp,
                                                    sensors, bitsize, 32, sim_repeats)

        # Run all distance sims
        if do_distance:
            print("Running simulations for distance/layout plot...")
            for layout in layouts:
                in_fpb = "input/layout_" + layout + "_sim_%03d_sensor%s.txt"
                out_fpb = "output/layout_" + layout + "_nav_%03d.txt"
                out_times_fp = "output/layout_" + layout + "_nav_times.txt"
                out_eif_fpb = "output/eif_layout_" + layout + "_nav_%03d.txt"

                if not DISTANCE_ONLY_EIF_BASE:
                    runner.run_simulation_repeats("input/track1.txt",
                                                    in_fpb,
                                                    out_fpb,
                                                    out_times_fp,
                                                    4, 1024, 32, sim_repeats)

                base_runner.run_normal_eif("input/track1.txt", in_fpb, out_eif_fpb, sim_repeats, 4)

        print("Finished running simulations.\n")


    """
     .d8888b.
    d88P  Y88b
    888    888
    888         .d88b.  88888b.d88b.  88888b.            .d88b.  888d888 888d888 .d88b.  888d888 .d8888b
    888        d88""88b 888 "888 "88b 888 "88b          d8P  Y8b 888P"   888P"  d88""88b 888P"   88K
    888    888 888  888 888  888  888 888  888          88888888 888     888    888  888 888     "Y8888b.
    Y88b  d88P Y88..88P 888  888  888 888 d88P d8b      Y8b.     888     888    Y88..88P 888          X88
     "Y8888P"   "Y88P"  888  888  888 88888P"  Y8P       "Y8888  888     888     "Y88P"  888      88888P'
                                      888
                                      888
                                      888
    """

    if do_sim_evaluate:
        # Compute encoding errors
        if do_encoding:
            print("Computing errors for encoding plot...")
            for frac_bits in encodings_to_test:
                out_fpb = "output/encoding_" + str(frac_bits) + "_nav_%03d.txt"
                errout_fpb = "output_evaluation/encoding_" + str(frac_bits) + "_nav_errors_%03d.txt"
                meanerrout_fp = "output_evaluation/encoding_" + str(frac_bits) + "_nav_mean_errors.txt"

                if not ENCODING_ONLY_EIF_BASE:
                    evaluator.create_sim_error_files("input/track1.txt", out_fpb, errout_fpb, meanerrout_fp, sim_repeats)
            
            # Only one EIF was run on the first encoding settings, as they don't affect the normal filter
            frac_bits = encodings_to_test[0]
            out_eif_fpb = "output/eif_encoding_" + str(frac_bits) + "_nav_%03d.txt"
            errout_eif_fpb = "output_evaluation/eif_encoding_" + str(frac_bits) + "_nav_errors_%03d.txt"
            meanerrout_eif_fp = "output_evaluation/eif_encoding_" + str(frac_bits) + "_nav_mean_errors.txt"
            evaluator.create_sim_error_files("input/track1.txt", out_eif_fpb, errout_eif_fpb, meanerrout_eif_fp, sim_repeats)

        # Compute timing errors
        if do_timing:
            print("Computing errors for timing plot...")
            for bitsize in timing_paillier_bitsizes_to_test:
                for sensors in range(min_timing_sensors_to_test, max_timing_sensors_to_test+1):
                    out_fpb = "output/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_%03d.txt"
                    errout_fpb = "output_evaluation/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_errors_%03d.txt"
                    meanerrout_fp = "output_evaluation/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_mean_errors.txt"
                    evaluator.create_sim_error_files("input/track1.txt", out_fpb, errout_fpb, meanerrout_fp, sim_repeats)

        # Compute distance errors
        if do_distance:
            print("Computing errors for distance/layout plot...")
            for layout in layouts:
                out_fpb = "output/layout_" + layout + "_nav_%03d.txt"
                errout_fpb = "output_evaluation/layout_" + layout + "_nav_errors_%03d.txt"
                meanerrout_fp = "output_evaluation/layout_" + layout + "_nav_mean_errors.txt"

                out_eif_fpb = "output/eif_layout_" + layout + "_nav_%03d.txt"
                errout_eif_fpb = "output_evaluation/eif_layout_" + layout + "_nav_errors_%03d.txt"
                meanerrout_eif_fp = "output_evaluation/eif_layout_" + layout + "_nav_mean_errors.txt"

                if not DISTANCE_ONLY_EIF_BASE:
                    evaluator.create_sim_error_files("input/track1.txt", out_fpb, errout_fpb, meanerrout_fp, sim_repeats)
                evaluator.create_sim_error_files("input/track1.txt", out_eif_fpb, errout_eif_fpb, meanerrout_eif_fp, sim_repeats)

        print("Finished computing errors.\n")

    """
    8888888b.  888          888
    888   Y88b 888          888
    888    888 888          888
    888   d88P 888  .d88b.  888888
    8888888P"  888 d88""88b 888
    888        888 888  888 888
    888        888 Y88..88P Y88b.
    888        888  "Y88P"   "Y888



    """

    if do_plot_create:
        # Initialise matplotlib for plotting or showing
        plotter.init_matplotlib_params(SAVE_NOT_SHOW_FIG, SHOW_LATEX_FIG)

        # Make encoding plot
        if do_encoding:
            print("Making encoding plot...")
            if encoding_plot_eif_base:
                plot_eif = 'output_evaluation/eif_encoding_%d_nav_mean_errors.txt' % encodings_to_test[0]
            else:
                plot_eif = None

            plotter.create_encoding_plots('output_evaluation/encoding_%d_nav_mean_errors.txt', encodings_to_test, 50, encoding_fig_width, plot_eif)

        # making timing plot
        if do_timing:
            print("Making timing plot...")
            plotter.create_timing_plots('output/timing_%d_%d_nav_times.txt', list(range(min_timing_sensors_to_test, max_timing_sensors_to_test+1)), timing_paillier_bitsizes_to_test, timing_fig_width)

        # Make layout distance plot
        if do_distance:
            print("Making distance/layout plot...")
            if layouts_plot_eif_base:
                plot_eif = 'output_evaluation/eif_layout_%s_nav_mean_errors.txt'
            else:
                plot_eif = None

            plotter.create_distance_plots('output_evaluation/layout_%s_nav_mean_errors.txt', layouts, layout_labels, 50, layouts_fig_width, plot_eif)

        print("Finished making plots.\n")

    return


if __name__ == '__main__':
    run_all(SIM_REPEATS, DO_MEASUREMENT_GEN, DO_SIM_RUN, DO_SIM_EVALUATE, DO_PLOT_CREATE, DO_ENCODING, DO_TIMING, DO_DISTANCE)