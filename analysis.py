import math
from functools import cached_property

import numpy as np


class SAvg:
    def __init__(self, s_up: float, s_down: float, s_up_uncertainty: float, s_down_uncertainty: float):
        self._s_up = s_up
        self._s_down = s_down
        self._s_up_uncertainty = s_up_uncertainty
        self._s_down_uncertainty = s_down_uncertainty

    @property
    def s_avg(self):
        weighted_sum = self._s_up / self._s_up_uncertainty ** 2 + self._s_down / self._s_down_uncertainty ** 2
        normalization = 1 / self._s_up_uncertainty ** 2 + 1 / self._s_down_uncertainty ** 2
        return weighted_sum / normalization

    @property
    def s_avg_delta(self):
        normalization = 1 / self._s_up_uncertainty ** 2 + 1 / self._s_down_uncertainty ** 2
        return math.sqrt(1 / normalization)


class ECalculator:
    def __init__(self, q1: float, q2: float, q1_delta: float, q2_delta: float):
        self._q1 = q1
        self._q2 = q2
        self._q1_delta = q1_delta
        self._q2_delta = q2_delta

    @cached_property
    def n(self):
        start = 4
        end = 12
        length = end - start + 1
        return np.linspace(start, end, length)

    @property
    def q1_per_n(self):
        return self._q1 / self.n

    @property
    def q1_per_n_delta(self):
        return self._q1_delta / self.n

    @property
    def q2_per_n(self):
        return self._q2 / self.n

    @property
    def q2_per_n_delta(self):
        return self._q2_delta / self.n
