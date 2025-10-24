import math

from experiment_calculations import ExperimentCalculations
from data import droplet_1_v0_times, droplet_1_v_up_times, droplet_1_v_down_times, droplet_1_eta, droplet_1_delta_eta
from data import droplet_2_v0_times, droplet_2_v_up_times, droplet_2_v_down_times, droplet_2_eta, droplet_2_delta_eta
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from util import chi_squared, linear_fit, sci
from analysis import SAvg, ECalculator
import numpy as np


droplet_1_experiment = ExperimentCalculations(droplet_1_v0_times, droplet_1_v_down_times, droplet_1_v_up_times, droplet_1_eta, droplet_1_delta_eta)
droplet_2_experiment = ExperimentCalculations(droplet_2_v0_times, droplet_2_v_down_times, droplet_2_v_up_times, droplet_2_eta, droplet_2_delta_eta)

v_terminal_1, v_terminal_un_1, v_terminal_exp_1 = sci(droplet_1_experiment.v_terminal, droplet_1_experiment.v_terminal_delta)
v_terminal_2, v_terminal_un_2, v_terminal_exp_2 = sci(droplet_2_experiment.v_terminal, droplet_2_experiment.v_terminal_delta)
print(rf"Droplet 1: $v_{{terminal}}$ = {v_terminal_1} $\pm$ {v_terminal_un_1} $m/s * 10^{v_terminal_exp_1:.0f}$")
print(rf"Droplet 2: $v_{{terminal}}$ = {v_terminal_2:.0f} $\pm$ {v_terminal_un_2:.0f} $m/s * 10^{v_terminal_exp_2:.0f}$")

dr_1, dr_un_1, dr_exp_1 = sci(droplet_1_experiment.drop_radius, droplet_1_experiment.drop_radius_delta)
dr_2, dr_un_2, dr_exp_2 = sci(droplet_2_experiment.drop_radius, droplet_2_experiment.drop_radius_delta)
print(rf"Droplet 1: $r_{{drop}}$ = {dr_1} $\pm$ {dr_un_1} $m * 10^{dr_exp_1:.0f}$")
print(rf"Droplet 2: $r_{{drop}}$ = {dr_2} $\pm$ {dr_un_2} $m * 10^{dr_exp_2:.0f}$")

def apply_curve_fit(experiment: ExperimentCalculations):
    popt_up, pcov_up = curve_fit(linear_fit, experiment.v_up_e_field, experiment.v_up,
                           sigma=experiment.v_up_delta, p0=[1, experiment.v_terminal])
    print(f"a = {popt_up[0]} +/- {np.sqrt(pcov_up[0][0])}, b = {popt_up[1]} +/- {np.sqrt(pcov_up[1][1])}")
    popt_down, pcov_down = curve_fit(linear_fit, experiment.v_down_e_field, experiment.v_down,
                           sigma=experiment.v_down_delta, p0=[1, experiment.v_terminal])
    print(f"a = {popt_down[0]} +/- {np.sqrt(pcov_down[0][0])}, b = {popt_down[1]} +/- {np.sqrt(pcov_down[1][1])}")
    return popt_up, pcov_up, popt_down, pcov_down

def set_axis_titles(ax: plt.Axes):
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x * 1e-3:.0f}'))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y * 1e4:.2f}'))
    ax.set_xlabel(r'Driving E-field Intensity ($\frac{kV}{m}$)')
    ax.set_ylabel(r'Velocity ($\frac{\mu m}{s}$)')

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

figsize=(8, 6)
margins = {
    'left': 0.12,
    'right': 0.89,
    'top': 0.87,
    'bottom': 0.12
}

# plt.figure(figsize=figsize)
fig, ax = plt.subplots(1, 1, figsize=figsize)
plt.subplots_adjust(**margins)

ax.errorbar(droplet_1_experiment.v_up_e_field, droplet_1_experiment.v_up, xerr=droplet_1_experiment.v_up_e_field_delta, yerr=droplet_1_experiment.v_up_delta, ecolor='red', capsize=2, fmt='x', label='v_up')
ax.errorbar(droplet_1_experiment.v_down_e_field, droplet_1_experiment.v_down, xerr=droplet_1_experiment.v_down_e_field_delta, yerr=droplet_1_experiment.v_down_delta, ecolor='red', capsize=2, fmt='x', label='v_down')
ax.plot(droplet_1_experiment.e_field_axis, linear_fit(droplet_1_experiment.e_field_axis, *popt_up_1), label=rf'v_up fit: $\chi^2/ndf = {up_1_chi_sq:.1f}/{up_1_ndf}$')
ax.plot(droplet_1_experiment.e_field_axis, linear_fit(droplet_1_experiment.e_field_axis, *popt_down_1), label=rf'v_down fit: $\chi^2/ndf = {down_1_chi_sq:.1f}/{down_1_ndf}$')
set_axis_titles(ax)

s_up_1, s_up_un_1, s_power_up_1 = sci(popt_up_1[0], np.sqrt(pcov_up_1[0][0]))
v0_up_1, v0_up_un_1, v0_power_up_1 = sci(popt_up_1[1], np.sqrt(pcov_up_1[1][1]))
s_down_1, s_down_un_1, s_power_down_1 = sci(popt_down_1[0], np.sqrt(pcov_down_1[0][0]))
v0_down_1, v0_down_un_1, v0_power_down_1 = sci(popt_down_1[1], np.sqrt(pcov_down_1[1][1]))

print("Drop 1:")
print(rf"s_up: {s_up_1} $\pm$ {s_up_un_1} $\frac{{m^2}}{{sV}} * 10^{{{s_power_up_1:.0f}}}$")
print(rf"v0_up: {v0_up_1} $\pm$ {v0_up_un_1} $m/s * 10^{{{v0_power_up_1:.0f}}}$")
print(rf"s_down: {s_down_1} $\pm$ {s_down_un_1} $\frac{{m^2}}{{sV}} * 10^{{{s_power_down_1:.0f}}}$")
print(rf"v0_down: {v0_down_1} $\pm$ {v0_down_un_1} $m/s * 10^{{{v0_power_down_1:.0f}}}$")

ax.set_title("Droplet 1 Velocity vs. E-field Intensity\n"
             rf"$s_{{v_{{up}}}}$ = {s_up_1} $\pm$ {s_up_un_1} $\frac{{m^2}}{{sV}}10^{{{s_power_up_1:.0f}}}$, "
             rf"$v_0$ = {v0_up_1} $\pm$ {v0_up_un_1} $\frac{{m}}{{s}}10^{{{v0_power_up_1:.0f}}}$;"
             rf"$s_{{v_{{down}}}}$ = {s_down_1} $\pm$ {s_down_un_1} $\frac{{m^2}}{{sV}}10^{{{s_power_down_1:.0f}}}$, "
             rf"$v_0$ = {v0_down_1} $\pm$ {v0_down_un_1} $\frac{{m}}{{s}}10^{{{v0_power_down_1:.0f}}}$")
plt.legend()
plt.show()

plt.close()

# plt.figure(figsize=figsize)
fig, ax = plt.subplots(1, 1, figsize=figsize)
plt.subplots_adjust(**margins)

ax.errorbar(droplet_2_experiment.v_up_e_field, droplet_2_experiment.v_up, xerr=droplet_2_experiment.v_up_e_field_delta, yerr=droplet_2_experiment.v_up_delta, ecolor='red', capsize=2, fmt='x', label='v_up')
ax.errorbar(droplet_2_experiment.v_down_e_field, droplet_2_experiment.v_down, xerr=droplet_2_experiment.v_down_e_field_delta, yerr=droplet_2_experiment.v_down_delta, ecolor='red', capsize=2, fmt='x', label='v_down')
ax.plot(droplet_2_experiment.e_field_axis, linear_fit(droplet_2_experiment.e_field_axis, *popt_up_2), label=rf'v_up fit: $\chi^2/ndf = {up_2_chi_sq:.1f}/{up_2_ndf}$')
ax.plot(droplet_2_experiment.e_field_axis, linear_fit(droplet_2_experiment.e_field_axis, *popt_down_2), label=rf'v_down fit: $\chi^2/ndf = {down_2_chi_sq:.1f}/{down_2_ndf}$')
set_axis_titles(ax)

s_up_2, s_up_un_2, s_power_up_2 = sci(popt_up_2[0], np.sqrt(pcov_up_2[0][0]))
v0_up_2, v0_up_un_2, v0_power_up_2 = sci(popt_up_2[1], np.sqrt(pcov_up_2[1][1]))
s_down_2, s_down_un_2, s_power_down_2 = sci(popt_down_2[0], np.sqrt(pcov_down_2[0][0]))
v0_down_2, v0_down_un_2, v0_power_down_2 = sci(popt_down_2[1], np.sqrt(pcov_down_2[1][1]))

print("Drop 2:")
print(rf"s_up: {s_up_2} $\pm$ {s_up_un_2} $\frac{{m^2}}{{sV}} * 10^{{{s_power_up_2:.0f}}}$")
print(rf"v0_up: {v0_up_2} $\pm$ {v0_up_un_2} $m/s * 10^{{{v0_power_up_2:.0f}}}$")
print(rf"s_down: {s_down_2} $\pm$ {s_down_un_2} $\frac{{m^2}}{{sV}} * 10^{{{s_power_down_2:.0f}}}$")
print(rf"v0_down: {v0_down_2} $\pm$ {v0_down_un_2} $m/s * 10^{{{v0_power_down_2:.0f}}}$")

ax.set_title("Droplet 2 Velocity vs. E-field Intensity\n"
             rf"$s_{{v_{{up}}}}$ = {s_up_2} $\pm$ {s_up_un_2} $\frac{{m^2}}{{sV}}10^{{{s_power_up_2:.0f}}}$, "
             rf"$v_0$ = {v0_up_2} $\pm$ {v0_up_un_2} $\frac{{m}}{{s}}10^{{{v0_power_up_2:.0f}}}$;"
             rf"$s_{{v_{{down}}}}$ = {s_down_2} $\pm$ {s_down_un_2} $\frac{{m^2}}{{sV}}10^{{{s_power_down_2:.0f}}}$, "
             rf"$v_0$ = {v0_down_2} $\pm$ {v0_down_un_2} $\frac{{m}}{{s}}10^{{{v0_power_down_2:.0f}}}$")

plt.legend()
plt.show()

# Slope calculations

s_1_avg_calc = SAvg(popt_up_1[0], popt_down_1[0], math.sqrt(pcov_up_1[0][0]), math.sqrt(pcov_down_1[0][0]))
s_2_avg_calc = SAvg(popt_up_2[0], popt_down_2[0], math.sqrt(pcov_up_2[0][0]), math.sqrt(pcov_down_2[0][0]))

s_1_avg, s_1_avg_un, s_1_power = sci(s_1_avg_calc.s_avg, s_1_avg_calc.s_avg_delta)
s_2_avg, s_2_avg_un, s_2_power = sci(s_2_avg_calc.s_avg, s_2_avg_calc.s_avg_delta)

print(rf"Droplet 1: $s_{{avg}}$ = {s_1_avg} $\pm$ {s_1_avg_un} $\frac{{m^2}}{{sV}} * 10^{{{s_1_power:.0f}}}$")
print(rf"Droplet 2: $s_{{avg}}$ = {s_2_avg} $\pm$ {s_2_avg_un} $\frac{{m^2}}{{sV}} * 10^{{{s_2_power:.0f}}}$")

q1 = droplet_1_experiment.calculate_q(s_1_avg_calc)
q1_delta = droplet_1_experiment.calculate_q_delta(s_1_avg_calc)

q2 = droplet_2_experiment.calculate_q(s_2_avg_calc)
q2_delta = droplet_2_experiment.calculate_q_delta(s_2_avg_calc)

q1_d, q1_unc_d, q1_exp = sci(q1, q1_delta)
q2_d, q2_unc_d, q2_exp = sci(q2, q2_delta)

print(rf"Droplet 1: $q_1$ = {q1_d} $\pm$ {q1_unc_d} C $10^{{{q1_exp}}}$")
print(rf"Droplet 2: $q_2$ = {q2_d}$ $\pm$ {q2_unc_d} C $10^{{{q2_exp}}}$")

print(rf"Droplet 1: $\frac{{q_1}}{{n}}$ = {q1_d} $\pm$ {q1_unc_d} C $10^{{{q1_exp}}}$")
print(rf"Droplet 2: $\frac{{q_2}}{{n}}$ = {q2_d}$ $\pm$ {q2_unc_d} C $10^{{{q2_exp}}}$")

plt.close()
fig, ax = plt.subplots()
e_calculator = ECalculator(q1, q2, q1_delta, q2_delta)
ax.errorbar(e_calculator.n, e_calculator.q1_per_n, yerr=e_calculator.q1_per_n_delta, fmt='x', capsize=7, ecolor='orange', elinewidth=0.7, label=r"$\frac{q_1}{n}$")
ax.errorbar(e_calculator.n, e_calculator.q2_per_n, yerr=e_calculator.q2_per_n_delta, fmt='x', capsize=7, ecolor='red', elinewidth=0.7, label=r"$\frac{q_2}{n}$")
ax.set_title("Integer Fractions of Droplet Charges")
ax.set_xlabel(r"Integer Divisor $n$")
ax.set_ylabel(r"Integer Fraction of Charge $\frac{q}{n}$")
plt.legend()
plt.show()

e_value = SAvg(q1 / 5, q2 / 6, q1_delta / 5, q2_delta / 6)
e_d, e_unc_d, e_exp = sci(e_value.s_avg, e_value.s_avg_delta)
print(rf"e- value: e: {e_d} $\pm$ {e_unc_d} C $10^{{{e_exp}}}$")

e_value = SAvg(q1 / 5, q2 / 7, q1_delta / 5, q2_delta / 7)
e_d, e_unc_d, e_exp = sci(e_value.s_avg, e_value.s_avg_delta)
print(rf"e- value: e: {e_d} $\pm$ {e_unc_d} C $10^{{{e_exp}}}$")

