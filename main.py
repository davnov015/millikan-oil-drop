from experiment_calculations import ExperimentCalculations
from data import droplet_1_v0_times, droplet_1_v_up_times, droplet_1_v_down_times, droplet_1_eta, droplet_1_delta_eta
from data import droplet_2_v0_times, droplet_2_v_up_times, droplet_2_v_down_times, droplet_2_eta, droplet_2_delta_eta
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from util import chi_squared, linear_fit
import numpy as np


droplet_1_experiment = ExperimentCalculations(droplet_1_v0_times, droplet_1_v_down_times, droplet_1_v_up_times, droplet_1_eta, droplet_1_delta_eta)
droplet_2_experiment = ExperimentCalculations(droplet_2_v0_times, droplet_2_v_down_times, droplet_2_v_up_times, droplet_2_eta, droplet_2_delta_eta)

def apply_curve_fit(experiment: ExperimentCalculations):
    popt_up, pcov_up = curve_fit(linear_fit, experiment.v_up_e_field, experiment.v_up,
                           sigma=experiment.v_up_delta, p0=[1, experiment.v_terminal])
    print(f"a = {popt_up[0]} +/- {np.sqrt(pcov_up[0][0])}, b = {popt_up[1]} +/- {np.sqrt(pcov_up[1][1])}")
    popt_down, pcov_down = curve_fit(linear_fit, experiment.v_down_e_field, experiment.v_down,
                           sigma=experiment.v_down_delta, p0=[1, experiment.v_terminal])
    print(f"a = {popt_down[0]} +/- {np.sqrt(pcov_down[0][0])}, b = {popt_down[1]} +/- {np.sqrt(pcov_down[1][1])}")
    return popt_up, pcov_up, popt_down, pcov_down

def set_axis_titles(ax: plt.Axes):
    ax.set_xlabel('Driving E-field (V/m)')
    ax.set_ylabel('Velocity (m/s)')

popt_up_1, pcov_up_1, popt_down_1, pcov_down_1 = apply_curve_fit(droplet_1_experiment)
popt_up_2, pcov_up_2, popt_down_2, pcov_down_2 = apply_curve_fit(droplet_2_experiment)

up_1_chi_sq = chi_squared(droplet_1_experiment.v_up, linear_fit(droplet_1_experiment.v_up_e_field, *popt_up_1), droplet_1_experiment.v_up_delta)
up_2_chi_sq = chi_squared(droplet_2_experiment.v_up, linear_fit(droplet_2_experiment.v_up_e_field, *popt_up_2), droplet_2_experiment.v_up_delta)
down_1_chi_sq = chi_squared(droplet_1_experiment.v_down, linear_fit(droplet_1_experiment.v_down_e_field, *popt_down_1), droplet_1_experiment.v_down_delta)
down_2_chi_sq = chi_squared(droplet_2_experiment.v_down, linear_fit(droplet_2_experiment.v_down_e_field, *popt_down_2), droplet_2_experiment.v_down_delta)

up_1_ndf = len(droplet_1_experiment.v_up) - 2
up_2_ndf = len(droplet_2_experiment.v_up) - 2
down_1_ndf = len(droplet_1_experiment.v_down) - 2
down_2_ndf = len(droplet_2_experiment.v_down) - 2

fig, ax = plt.subplots(1, 1)
plt.subplots_adjust(left=0.15, right=0.97, top=0.92, bottom=0.12)

ax.errorbar(droplet_1_experiment.v_up_e_field, droplet_1_experiment.v_up, xerr=droplet_1_experiment.v_up_e_field_delta, yerr=droplet_1_experiment.v_up_delta, fmt='x', label='v_up')
ax.errorbar(droplet_1_experiment.v_down_e_field, droplet_1_experiment.v_down, xerr=droplet_1_experiment.v_down_e_field_delta, yerr=droplet_1_experiment.v_down_delta, fmt='x', label='v_down')
ax.plot(droplet_1_experiment.e_field_axis, linear_fit(droplet_1_experiment.e_field_axis, *popt_up_1), label=rf'v_up fit: $\chi^2/ndf = {up_1_chi_sq:.2f}/{up_1_ndf}$')
ax.plot(droplet_1_experiment.e_field_axis, linear_fit(droplet_1_experiment.e_field_axis, *popt_down_1), label=rf'v_down fit: $\chi^2/ndf = {down_1_chi_sq:.2f}/{down_1_ndf}$')
set_axis_titles(ax)
ax.set_title("Droplet 1 Velocity vs. E-field")
plt.legend()
plt.show()

plt.close()

fig, ax = plt.subplots(1, 1)
plt.subplots_adjust(left=0.15, right=0.97, top=0.92, bottom=0.12)

ax.errorbar(droplet_2_experiment.v_up_e_field, droplet_2_experiment.v_up, xerr=droplet_2_experiment.v_up_e_field_delta, yerr=droplet_2_experiment.v_up_delta, fmt='x', label='v_up')
ax.errorbar(droplet_2_experiment.v_down_e_field, droplet_2_experiment.v_down, xerr=droplet_2_experiment.v_down_e_field_delta, yerr=droplet_2_experiment.v_down_delta, fmt='x', label='v_down')
ax.plot(droplet_2_experiment.e_field_axis, linear_fit(droplet_2_experiment.e_field_axis, *popt_up_2), label=rf'v_up fit: $\chi^2/ndf = {up_2_chi_sq:.2f}/{up_2_ndf}$')
ax.plot(droplet_2_experiment.e_field_axis, linear_fit(droplet_2_experiment.e_field_axis, *popt_down_2), label=rf'v_down fit: $\chi^2/ndf = {down_2_chi_sq:.2f}/{down_2_ndf}$')
set_axis_titles(ax)
ax.set_title("Droplet 2 Velocity vs. E-field")
plt.legend()
plt.show()
