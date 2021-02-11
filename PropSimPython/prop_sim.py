# PropSim in Python 3.8
# Project Caelus, Aphlex 1C Engine
# 09 February, 2021
# Authors: Jason Chen, Liam West, Anya Mischel, Jack Blair

import numpy as np
import time
from scipy.integrate import solve_ivp
from typing import Tuple
from helpers.n2o import n2o_properties, n2o_find_T
from helpers.state_flow import liquid_state_vec
from helpers.wrappers import fast_interp_1, fast_interp_2, fast_interp_3
from main_functs import integration, find_G, find_mass_gradient


#----------Unit Conversions-------------
psi_to_Pa = 6894.75729 # 1 psi in Pa
in_to_m = 0.0254 # 1 in in m
mm_to_m = 1e-3 # 1 mm in m
lbf_to_N = 4.44822162 # 1 lbf in N
lbm_to_kg = 0.453592 # 1 lbm in kg
atm_to_Pa = 101325 # 1 atm in Pa
L_to_m3 = 1e-3 # 1 L in m^3

inputs = {
    #-------------Constants-------------
    "constants": {
        "R_u": 8.3144621, # Universal gas constant [J/mol*K]
        "g_0": 9.80665 # Standard gravitational constant [m/s^2]
    },

    #---------Unit Conversions----------
    "unit_conv": {
        "psi_to_Pa": 6894.75729,
        "in_to_m": 0.0254,
        "mm_to_m": 1e-3,
        "lbf_to_N": 4.44822162,
        "lbm_to_kg": 0.453592,
        "atm_to_Pa": 101325,
        "L_to_m3": 1e-3
    },

    #--------Runtime Properties---------
    "options": {
        # Simulation time limit
        "t_final": 60, # sec
        # Output time resolution
        "dt": 0.01, # sec
        # True for plots, False for no plots
        "output_on": True
    },

    #--------Oxidizer Properties---------
    "ox": {
        "injector_area": 2.571e-05, # m^2
        # Discharge coefficient (Cd)
        "injector_Cd": 0.9, 
        # Tank volume
        "tank_V": 11.564*L_to_m3, # m^3
        # Tank temperature
        "tank_T": 300, # K
        # Tank inner diameter (ID)
        "tank_id": 5*in_to_m,
        # Distance from bottom of tank to injector
        "tank_h_offset": 2, # m
        # Nitrous volume
        "liquid_V": 4.208*L_to_m3,
        # Main flow line diameter
        "d_flowline": .5*in_to_m,
    },

    #----------Fuel Properties-----------
    "fuel": {
        "injector_area": 6.545e-06, # m^2
        # Discharge coefficient (Cd)
        "injector_Cd": 0.88, 
        # Tank volume
        "tank_V": 11.564*L_to_m3, # m^3
        # Tank temperature
        "tank_T": 300, # K
        # Tank inner diameter (ID)
        "tank_id": 5*in_to_m,
        # Distance from bottom of tank to injector
        "tank_h_offset": 24*in_to_m, # m
        # Nitrous volume
        "liquid_V": 1.267*L_to_m3,
        # Main flow line diameter
        "d_flowline": .5*in_to_m,
        # Propellant density
        "rho": 789
    },

    #-----Fuel Pressurant Properties------
    "fuel_pressurant": {
        "gas_props": {
            # Specific volume
            "c_v": None,
            # Molecular mass
            "mol_mass": None
        },
        # The items below aren't used in pressure blowdown
        "set_pressure": 350*psi_to_Pa,
        "storage_init_press": 4500*psi_to_Pa,
        "stroage_tank_V": 0.0*L_to_m3,
        "flow_CdA": 8*mm_to_m**2 # m^2
    },

    #----------Other Properties-----------
    # Ball valve time to injector area (s)
    "dt_valve_open": 0.01, # sec
    # Rocket dry mass
    "rocket_dry_mass": 30*lbm_to_kg, # kg
    # Combustion chamber dimensions
    "length_cc": 0.258, # m
    "d_cc": 0.08, # m
    # Estimated nozzle efficiency
    "nozzle_eff": 0.95,
    # Nozzle correction factor
    "nozzle_correction_factor": 0.9830,
    # Estimated combustion efficiency
    "c_star_efficiency": 0.85,
    # Nozzle throat diameter
    "d_throat": 30.46e-3,
    # Expansion Ratio
    "exp_ratio": 3.26,
    # Ambient Temperature
    "T_amb": 292, # K
    # Ambient Pressure
    "p_amb": 9.554e04, # Pa
}


def run_performance():
    """
    Runs integration.py to integrate the state vector and records output
    """
    # Get oxidizer properties at the given temperature
    n2o = n2o_properties(inputs["ox"]["T_tank"])
    # Our integration variables are oxidizer mass and liquid oxidizer volume
    Mox = n2o["rho_l"]*(inputs["ox"]["liquid_V"]) + n2o["rho_g"]*(inputs["ox"]["tank_V"]-inputs["ox"]["liquid_V"])
    if inputs["options"]["output_on"]:
        print("Initial oxidizer mass: {} kg.".format(Mox))
    start = time.perf_counter() # Start timer for integration

    time, record = integration(inputs) # Time = time for integration, record = output data
    F_thrust        = record.F_thrust
    p_cc            = record.p_cc
    p_oxtank        = record.p_oxtank
    p_oxpresstank   = record.p_oxpresstank
    p_fueltank      = record.p_fueltank
    p_fuelpresstank = record.p_fuelpresstank
    p_oxmanifold    = record.p_oxmanifold
    T_oxtank        = record.T_oxtank
    T_cc            = record.T_cc
    area_core       = record.area_core
    OF              = record.OF_i
    gamma_ex        = record.gamma_ex
    m_dot_ox        = record.m_dot_ox
    m_dot_fuel      = record.m_dot_fuel
    p_crit          = record.p_crit
    m_dot_ox_crit   = record.m_dot_ox_crit
    M_e             = record.M_e
    p_exit          = record.p_exit
    p_shock         = record.p_shock

    time_elapsed = start-time.perf_counter()
    if inputs["options"]["output_on"]:
        print("Time elapsed for this timestep: {} sec.".format(time_elapsed))


if __name__ == "__main__":
    run_performance()
