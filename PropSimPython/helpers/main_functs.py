# Helper Methods for PropSimPython
# Project Caelus, Aphlex 1C Engine
# 02 February, 2021

import numpy as np
from typing import Tuple


def integration(inputs: dict) -> Tuple[float, dict]:
    record = {
        "F_thrust": None,
        "p_cc": None,
        "p_oxtank": None,
        "p_oxpresstank": None,
        "p_fueltank": None,
        "p_fuelpresstank": None,
        "p_oxmanifold": None,
        "T_oxtank": None,
        "T_cc": None,
        "area_core": None,
        "OF_i": None,
        "gamma_ex": None,
        "m_dot_ox": None,
        "m_dot_fuel": None,
        "p_crit": None,
        "m_dot_ox_crit": None,
        "M_e": None,
        "p_exit": None,
        "p_shock": None,
    } # Recording the output data for this timestep

    state_0, x0 = InitializeLiquidState(inputs)
    pass


def find_G():
    raise NotImplementedError


def find_mass_gradient():
    raise NotImplementedError
