from experiment_calculations import ExperimentCalculations
from data import droplet_1_v0_times, droplet_1_v_up_times, droplet_1_v_down_times, droplet_1_eta, droplet_1_delta_eta
from data import droplet_2_v0_times, droplet_2_v_up_times, droplet_2_v_down_times, droplet_2_eta, droplet_2_delta_eta
from matplotlib import pyplot as plt


droplet_1_experiment = ExperimentCalculations(droplet_1_v0_times, droplet_1_v_down_times, droplet_1_v_up_times, droplet_1_eta, droplet_1_delta_eta)
droplet_2_experiment = ExperimentCalculations(droplet_2_v0_times, droplet_2_v_down_times, droplet_2_v_up_times, droplet_2_eta, droplet_2_delta_eta)


fig, ax = plt.subplots(1, 1)
ax.errorbar(droplet_1_experiment.v_up_e_field, droplet_1_experiment.v_up, xerr=droplet_1_experiment.v_up_e_field_delta, yerr=droplet_1_experiment.v_up_delta, fmt='x', label='Upward velocity')
ax.errorbar(droplet_1_experiment.v_down_e_field, droplet_1_experiment.v_down, xerr=droplet_1_experiment.v_down_e_field_delta, yerr=droplet_1_experiment.v_down_delta, fmt='x', label='Downward velocity')
ax.set_xlabel('Driving E-field (V/m)')
ax.set_ylabel('Velocity (m/s)')
plt.legend()
plt.show()

plt.close()

fig, ax = plt.subplots(1, 1)
ax.errorbar(droplet_2_experiment.v_up_e_field, droplet_2_experiment.v_up, xerr=droplet_2_experiment.v_up_e_field_delta, yerr=droplet_2_experiment.v_up_delta, fmt='x', label='Upward velocity')
ax.errorbar(droplet_2_experiment.v_down_e_field, droplet_2_experiment.v_down, xerr=droplet_2_experiment.v_down_e_field_delta, yerr=droplet_2_experiment.v_down_delta, fmt='x', label='Downward velocity')
ax.set_xlabel('Driving E-field (V/m)')
ax.set_ylabel('Velocity (m/s)')
plt.legend()
plt.show()
