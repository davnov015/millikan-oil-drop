from functools import cached_property

from data import (
    vertical_displacement,
    delta_vertical_displacement,
    plate_separation,
    delta_plate_separation,
    delta_voltage
)
from analysis import SAvg
import numpy as np
import math

from typing import Union


class ExperimentCalculations:
    """
    Performs calculations on experimental data associated with a given trial.
    """

    vertical_displacement = vertical_displacement
    delta_vertical_displacement = delta_vertical_displacement
    plate_separation = plate_separation
    delta_plate_separation = delta_plate_separation
    delta_voltage = delta_voltage

    P = 101300  # barometric pressure, Pa
    B = 8.22e-3  # constant, Pa * m
    g = 9.81  # Unit: m/s^2 Gravitational acceleration.
    rho = 886  # Unit: kg/m^3 This is oil density.

    def __init__(
        self,
        v_terminal_times: list[float],
        v_down_times: dict[int, list[float]],
        v_up_times: dict[int, list[float]],
        eta: float,
        eta_delta: float,
    ):
        """

        :param v_terminal_times: The terminal velocity measurement times.
        :param v_down_times: The downward velocity measurement times, with driving E-field.
        :param v_up_times: The upward velocity measurement times, with driving E-field.
        :param eta: The eta coefficient, based on chamber temperature (Units: Ns/m^2).
        :param eta_delta: The uncertainty in eta.
        """
        self._v_terminal_times = v_terminal_times
        self._v_down_times = v_down_times
        self._v_up_times = v_up_times
        self._eta = eta
        self._eta_delta = eta_delta

    @cached_property
    def v_terminal(self) -> float:
        """
        Calculates the signed terminal velocity. This is a negative quantity since the drop is falling.
        :return: The terminal velocity in m/s.
        """
        return -self._get_velocity(self._v_terminal_times)

    @cached_property
    def v_terminal_delta(self) -> float:
        """Calculates the uncertainty in the terminal velocity."""
        return self._get_velocity_delta(self._v_terminal_times)

    def __sqrt_term(self):
        return math.sqrt((self.B / (2 * self.P)) ** 2 - 9 * self._eta * self.v_terminal / (2 * self.g * self.rho))

    @cached_property
    def drop_radius(self) -> float:
        """
        Calculates the drop radius based on the terminal velocity.
        :return: The drop radius in meters.
        """
        assert self.v_terminal < 0, "v_terminal must be negative (i.e. the signed velocity of the falling drop)"
        return self.__sqrt_term() - self.B / (2 * self.P)

    @property
    def drop_mass(self):
        return 4 / 3 * math.pi * self.drop_radius ** 3 * self.rho

    @property
    def drop_mass_delta(self):
        dMdA = 4 * math.pi * self.drop_radius ** 2 * self.rho
        return dMdA * self.drop_radius_delta

    def calculate_q(self, s_avg: SAvg):
        return -s_avg.s_avg * self.drop_mass * self.g / self.v_terminal

    def calculate_q_delta(self, s_avg: SAvg):
        dQdS = -self.drop_mass * self.g / self.v_terminal
        dQdM = s_avg.s_avg * self.g / self.v_terminal
        dQdV0 = s_avg.s_avg * self.drop_mass * self.g / self.v_terminal
        return math.sqrt((dQdS * s_avg.s_avg_delta) ** 2 + (dQdM * self.drop_mass_delta) ** 2 + (dQdV0 * self.v_terminal_delta) ** 2)

    @cached_property
    def drop_radius_delta(self) -> float:
        """Calculates the uncertainty in the drop radius."""
        dRdV0 = 1 / (2 * self.__sqrt_term()) * (-9 * self._eta / (2 * self.g * self.rho))
        dRdEta = 1 / (2 * self.__sqrt_term()) * (-9 * self.v_terminal / (2 * self.g * self.rho))
        return math.sqrt((dRdV0 * self.v_terminal_delta) ** 2 + (dRdEta * self._eta_delta) ** 2)

    @property
    def v_up(self) -> np.ndarray:
        """
        Calculates the velocity of the rise when the driving E-field points upwards.
        :return: A list of velocities in m/s, for each E-field strength.
        """
        return self._v_up_calculation[:, 1]

    @property
    def v_up_delta(self) -> np.ndarray:
        """Calculates the uncertainty of the raise velocities when the driving E-field points upwards."""
        return self._v_up_calculation[:, 3]

    @property
    def v_up_e_field(self) -> np.ndarray:
        """Calculates the driving E-field strength when the drop is rising."""
        return self._v_up_calculation[:, 0]

    @property
    def v_up_e_field_delta(self) -> np.ndarray:
        """Calculates the uncertainty in the driving E-field strength when the drop is rising."""
        return self._v_up_calculation[:, 2]

    @property
    def v_down(self) -> np.ndarray:
        """
        Calculates the velocity of the drop when the driving E-field points downwards.
        :return: A list of velocities in m/s, for each E-field strength.
        """
        return self._v_down_calculation[:, 1]

    @property
    def v_down_delta(self) -> np.ndarray:
        """Calculates the uncertainty of the drop velocities when the driving E-field points downwards."""
        return self._v_down_calculation[:, 3]

    @property
    def v_down_e_field(self) -> np.ndarray:
        """Calculates the driving E-field strength when the drop is falling."""
        return self._v_down_calculation[:, 0]

    @cached_property
    def v_down_e_field_delta(self) -> np.ndarray:
        """Calculates the uncertainty in the driving E-field strength when the drop is falling."""
        return self._v_down_calculation[:, 2]

    @cached_property
    def e_field_axis(self):
        """Generate the E-field axis, suitable for plotting."""
        min_e_field = min(self.v_down_e_field.min(), self.v_up_e_field.min())
        max_e_field = max(self.v_down_e_field.max(), self.v_up_e_field.max())
        return np.linspace(min_e_field, max_e_field, 100)

    @cached_property
    def _v_up_calculation(self) -> np.ndarray:
        return self._v_calculation(self._v_up_times)

    @cached_property
    def _v_down_calculation(self) -> np.ndarray:
        return self._v_calculation(self._v_down_times)

    @classmethod
    def _v_calculation(cls, voltage_v_times_mapping: dict[int, list[float]]):
        """
        Calculates the velocities associated with different E-field strengths. The uncertainties are calculated.
        :param voltage_v_times_mapping: A mapping of voltage to velocity times.
        :return: A series of the velocity times and their uncertainties.
        """
        series = np.zeros([len(voltage_v_times_mapping), 4])
        for i, (voltage, v_times) in enumerate(voltage_v_times_mapping.items()):
            e_field = cls._e_field(voltage)
            series[i, 0] = e_field
            series[i, 1] = cls._get_velocity(v_times)
            series[i, 2] = cls._e_field_delta(voltage)
            series[i, 3] = cls._get_velocity_delta(v_times)
        series = series[series[:, 0].argsort()]
        assert np.all(np.diff(series[:, 0]) > 0), "v_times must be sorted by E-field values"
        return series

    @classmethod
    def _e_field(cls, voltage: float):
        """Calculates the E-field between the plates."""
        return voltage / cls.plate_separation

    @classmethod
    def _e_field_delta(cls, voltage: float):
        """Calculates the E-field uncertainty."""
        return (voltage / cls.plate_separation *
                math.sqrt((cls.delta_voltage / voltage) ** 2 + (cls.delta_plate_separation / cls.plate_separation) ** 2))

    @classmethod
    def _time_measurement_delta(cls, timings: Union[list[float], np.ndarray]):
        """Calculates the time measurement uncertainty. We use the standard deviation of the measurements."""
        return np.std(timings)

    @classmethod
    def _get_velocity(cls, timings: Union[list[float], np.ndarray]):
        """Calculates the velocity of the drop."""
        return cls.vertical_displacement / np.average(timings)

    @classmethod
    def _get_velocity_delta(cls, timings: Union[list[float], np.ndarray]):
        """Calculates the uncertainty in the velocity of the drop."""
        velocity = cls._get_velocity(timings)
        return velocity * math.sqrt((cls.delta_vertical_displacement / cls.vertical_displacement) ** 2 +
                                    (cls._time_measurement_delta(timings) / np.average(timings)) ** 2)
