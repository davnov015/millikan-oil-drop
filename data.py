import numpy as np

### Setup-related constants ###
# Note: variables prefixed with 'delta_' represent uncertainties

# The vertical displacement associated with velocity calculation (up/down times)
vertical_displacement = 0.001 # Unit: m
delta_vertical_displacement = 1 / 20 * vertical_displacement

plate_separation = 0.0076 # Unit: m
delta_plate_separation = 0.0001 # Unit: m

delta_voltage = 5   # Unit: V
delta_time = 0.500  # Unit s This is the total uncertainty of each time; but maybe use 1*sigma of the time distribution?

### Experimental data ###

# Droplet 1

droplet_1_eta = 1.859 * 1e-5 # Unit: Ns/m^2
droplet_1_delta_eta = 0.003 * 1e-5

droplet_1_v0_times = [13.15, 12.37, 10.40, 13.39]

droplet_1_v_down_times = {
    500: [3.21, 3.58, 3.18, 3.65],
    400: [3.73, 4.13, 3.63, 4.01],
    300: [4.26, 4.85, 4.81, 4.48],
    200: [5.23, 5.21, 5.66, 4.50],
    100: [7.85, 7.15, 7.68, 7.00],
}

droplet_1_v_up_times = {
    500: [7.23, 6.21, 7.42, 7.10],
    400: [10.98, 10.63, 10.90, 12.43],
    300: [18.25, 15.45, 16.00, 18.19],
    200: [24.51, 22.26, 27.13, 24.47],
}

# Droplet 1

droplet_2_eta = 1.862 * 1e-5   # Unit: Ns/m^2
droplet_2_delta_eta = 0.001 * 1e-5

droplet_2_v0_times = [11.29, 12.06, 12.07, 15.55]

droplet_2_v_up_times = {
    500: [7.20, 5.87, 4.91, 5.07],
    400: [5.99, 6.58, 5.92, 5.53],
    300: [7.45, 7.97, 7.77, 9.70],
    200: [16.77, 20.96, 20.77, 21.10],
}

droplet_2_v_down_times = {
    500: [3.27, 3.07, 3.08, 2.75],
    400: [2.93, 3.05, 3.17, 2.67],
    300: [4.25, 3.80, 3.92, 4.23],
    200: [5.43, 4.25, 5.07, 5.55],
    100: [7.30, 6.43, 6.43, 6.17],
}
