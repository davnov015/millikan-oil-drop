import math


P = 101300  # barometric pressure, Pa
B = 8.22e-3 # constant, Pa * m
g = 9.81 # m/s^2
rho = 886 # kg/m^3
eta = 1.86e-5 # Ns/m^2

def get_drop_radius(v_terminal: float, eta: float) -> float:
    """
    Calculates the drop radius based on the terminal velocity.
    :param v_terminal: The signed terminal velocity - this must be a negative value.
    :param eta: The eta coefficient, based on chamber temperature (Units: Ns/m^2).
    :return: The drop radius in meters.
    """
    assert v_terminal < 0, "v_terminal must be negative (i.e. the signed velocity of the falling drop)"
    return math.sqrt((B / (2 * P)) ** 2 - (9 * eta * v_terminal / (2 * g * rho))) - B / (2 * P)
