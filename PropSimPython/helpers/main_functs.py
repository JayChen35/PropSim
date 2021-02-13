# Main helper methods for PropSimPython
# Project Caelus, Aphlex 1C Engine
# 02 February, 2021

import numpy as np
from typing import Tuple
from state_flow import init_liquid_state, LiquidStateVector


def integration(inputs: dict) -> Tuple[float, dict]:
    """
    Integrate necesary differential equations for rocket engine modeling using Euler's method. 
    TODO: Include real combustion properties and supercharging.

    INPUTS:
        - inputs: structure of motor characteristics (all units SI: m, s, 
        kg, K, mol)
            - CombustionData: string filename in "Combustion Data/" folder 
            from which to source combustion data
            - Injexit_area: total orifice area of injector
            - Cdischarge: discharge coefficient of injector orifices
            - MOVTime: time for injector flow to ramp up to 100
            - rocket_dry_mass: dry mass of rocket
            - tankvol: volume of oxidizer tank
            - l_vol: initial volume of liquid nitrous oxide in oxidizer 
            tank
            - tank_id: inner diameter of tank
            - h_offset: height difference between bottom of tank and 
            injector
            - flowline_id: inner diameter of flow line
            - Ttank: tank initial temperature
            - Fueldensity: density of fuel
            - grainlength: length of fuel grain
            - graindiameter: outer diameter of grain
            - portradius0: grain intial port radius
            - chamberlength: combustion chamber total length
            - n: grain ballistic coefficient (r_dot = a*G^n)
                - a: grain ballistic coefficient (r_dot = a*G^n)
            - nozzle_efficiency: nozzle exhaust energy efficiency 
            (proportion)
            - nozzle_correction_factor: proportion of ideal thrust 
            (factor for divergence losses, etc.) actually achieved
            - c_star_efficiency: combustion efficiency / proportion of c*
            actually achieved
            - Tdiameter: nozzle throat diameter
            - E: expansion ratio
            - Tamb: ambient temperature
            - Pamb: ambient pressure
            - SPress: supercharging regulator pressure
            - M_sc: molecular mass
            - SVol: volume of external pressurization tank (0 for 
            supercharging / no external tank)
            - c_v_S: specific heat at constant volume of pressurant gas
                - c_p_S: specific heat at constant volume of pressurant gas
            - P_S_i: initial pressurant gas storage pressure (must be 
            present, but not used for supercharging / no external tank)
            - S_CdA: effective flow area of pressurant gas (must be 
            but not used for supercharging / no external tank)
            - Charged: 1 for pressurant gas present, 0 for no pressurant 
            gas present
        - mode: structure of options defining mode of motor operation
            - combustion_on: 1 for hot-fire, 0 for cold-flow
            - flight_on: 1 for flight conditions, 0 for ground conditions
        - tspan: vector of time values over which to record outputs
    OUTPUS:
        - tspan: output time vector
        - F_thrust: thrust over tspan
        - p_cc: combustion chamber pressure over tspan
        - p_oxtank: tank pressure over tspan
        - p_oxmanifold: oxidizer manifold pressure over tspan
        - T_tank: tank temperature over tspan
        - T_cc: combustion chamber as temperature over tspan
        - area_core: fuel grain core area over tspan
        - OF: oxidizer/fuel ratio over tspan
        - m_dot_ox: oxidizer mass flow rate over tspan
        - m_dot_ox_crit: critical two-phase oxidizer mass flow rate over tspan
        - p_crit: two-phase critical downstream pressure over tspan
        - M_e: exit Mach number over tspan
        - p_exit: nozzle exit pressure over tspan
        - p_shock: critical back pressure for normal shock formation over
        tspan
    """
    # Recording the output data for this timestep
    record = { 
        "F_thrust"        : None,
        "p_cc"            : None,
        "p_oxtank"        : None,
        "p_oxpresstank"   : None,
        "p_fueltank"      : None,
        "p_fuelpresstank" : None,
        "p_oxmanifold"    : None,
        "T_oxtank"        : None,
        "T_cc"            : None,
        "area_core"       : None,
        "OF_i"            : None,
        "gamma_ex"        : None,
        "m_dot_ox"        : None,
        "m_dot_fuel"      : None,
        "p_crit"          : None,
        "m_dot_ox_crit"   : None,
        "M_e"             : None,
        "p_exit"          : None,
        "p_shock"         : None
    } 
    state_0, x0 = init_liquid_state(inputs)
    

def find_G():
    raise NotImplementedError


def find_mass_gradient():
    raise NotImplementedError
