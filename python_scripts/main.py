"""

"""

import generator
import runner
import evaluator
import plotter

# Process (gen., sim., eval., etc.) repeats 
SIM_REPEATS = 3

# Which part of process to do
DO_MEASUREMENT_GEN = True
DO_SIM_RUN = True
DO_SIM_EVALUATE = True
DO_PLOT_CREATE = True

# Which scenarios to run
DO_ENCODING = True
DO_TIMING = True
DO_DISTANCE = True

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

    # Encoding test params
    encodings_to_test = [(64, 16), (128, 32), (256, 64)]

    # Timing test params
    timing_paillier_bitsizes_to_test = [512, 1024, 2048]
    timing_sensors_to_test = 8

    # Distance test params
    layouts = ['normal', 'big', 'verybig', 'small']

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
            generator.generate_sim_inputs("input/track1.txt", "input/encoding_sim_%03d_sensor%s.txt", sim_repeats, 0, 4)

        # Genereate measurements for timing plots
        if do_timing:
            generator.generate_sim_inputs("input/track1.txt", "input/timing_sim_%03d_sensor%s.txt", sim_repeats, 0, 8)

        # Generate measurements for distance plots
        if do_distance:
            for i,layout in enumerate(layouts):
                sensor_fpb = "input/layout_" + layout + "_sim_%03d_sensor%s.txt"
                sensor_start_index = i*4
                generator.generate_sim_inputs("input/track1.txt", sensor_fpb, sim_repeats, sensor_start_index, 4)

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
            for mod_bits, frac_bits in encodings_to_test:

                out_fpb = "output/encoding_" + str(mod_bits) + "_" + str(frac_bits) + "_nav_%03d.txt"
                out_times_fp = "output/encoding_" + str(mod_bits) + "_" + str(frac_bits) + "_nav_times.txt"

                runner.run_simulation_repeats("input/track1.txt",
                                                "input/encoding_sim_%03d_sensor%s.txt",
                                                out_fpb,
                                                out_times_fp,
                                                4, 1024, mod_bits, frac_bits, sim_repeats)

        # Run all timing sims
        if do_timing:
            for bitsize in timing_paillier_bitsizes_to_test:
                for sensors in range(3, timing_sensors_to_test+1):

                    out_fpb = "output/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_%03d.txt"
                    out_times_fp = "output/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_times.txt"
                    runner.run_simulation_repeats("input/track1.txt",
                                                    "input/timing_sim_%03d_sensor%s.txt",
                                                    out_fpb,
                                                    out_times_fp,
                                                    sensors, bitsize, 128, 32, sim_repeats)

        # Run all distance sims
        if do_distance:
            for layout in layouts:
                in_fpb = "input/layout_" + layout + "_sim_%03d_sensor%s.txt"
                out_fpb = "output/layout_" + layout + "_nav_%03d.txt"
                out_times_fp = "output/layout_" + layout + "_nav_times.txt"
                runner.run_simulation_repeats("input/track1.txt",
                                                in_fpb,
                                                out_fpb,
                                                out_times_fp,
                                                4, 1024, 128, 32, sim_repeats)

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
            for mod_bits, frac_bits in encodings_to_test:
                out_fpb = "output/encoding_" + str(mod_bits) + "_" + str(frac_bits) + "_nav_%03d.txt"
                errout_fpb = "output_evaluation/encoding_" + str(mod_bits) + "_" + str(frac_bits) + "_nav_errors_%03d.txt"
                meanerrout_fp = "output_evaluation/encoding_" + str(mod_bits) + "_" + str(frac_bits) + "_nav_mean_errors.txt"
                evaluator.create_sim_error_files("input/track1.txt", out_fpb, errout_fpb, meanerrout_fp, sim_repeats)

        # Compute timing errors
        if do_timing:
            for bitsize in timing_paillier_bitsizes_to_test:
                for sensors in range(3, timing_sensors_to_test+1):
                    out_fpb = "output/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_%03d.txt"
                    errout_fpb = "output_evaluation/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_errors_%03d.txt"
                    meanerrout_fp = "output_evaluation/timing_" + str(bitsize) + "_" + str(sensors) + "_nav_mean_errors.txt"
                    evaluator.create_sim_error_files("input/track1.txt", out_fpb, errout_fpb, meanerrout_fp, sim_repeats)

        # Compute distance errors
        if do_timing:
            for layout in layouts:
                out_fpb = "output/layout_" + layout + "_nav_%03d.txt"
                errout_fpb = "output_evaluation/layout_" + layout + "_nav_errors_%03d.txt"
                meanerrout_fp = "output_evaluation/layout_" + layout + "_nav_mean_errors.txt"
                evaluator.create_sim_error_files("input/track1.txt", out_fpb, errout_fpb, meanerrout_fp, sim_repeats)

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

    print("TODO: plots")

    return


if __name__ == '__main__':
    run_all(SIM_REPEATS, DO_MEASUREMENT_GEN, DO_SIM_RUN, DO_SIM_EVALUATE, DO_PLOT_CREATE, DO_ENCODING, DO_TIMING, DO_DISTANCE)