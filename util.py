import numpy as np
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
