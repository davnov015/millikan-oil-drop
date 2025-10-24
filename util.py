import numpy as np
import math
from typing import Union


def chi_squared(data: np.ndarray, model: np.ndarray, sigma: np.ndarray) -> float:
    """
    Calculates the chi-squared value for a given data set.
    :param data: The data set.
    :param model: The modeled data set.
    :param sigma: The uncertainties on the data.
    :return: Chi-squared value.
    """
    return np.sum(((data - model) / sigma) ** 2)

def linear_fit(x: Union[float, np.ndarray], m: float, b: float) -> Union[float, np.ndarray]:
    """
    Linear fit function.
    :param x: Independent variable.
    :param m: Slope.
    :param b: Y-intercept.
    :return:
    """
    return m * x + b

def sci(number: float, uncertainty: float):
    assert number < 1 and uncertainty < 1, "Number and uncertainty must be less than 1"
    exponent = math.ceil(math.log10(math.fabs(1 / number)))
    number *= 10 ** exponent
    uncertainty *= 10 ** exponent
    unc_exponent = math.ceil(math.log10(math.fabs(1 / uncertainty)))
    uncertainty = round(uncertainty * 10 ** unc_exponent) / 10 ** unc_exponent
    number = round(number, unc_exponent)
    return number, uncertainty, -exponent
