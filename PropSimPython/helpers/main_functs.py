# Main helper methods for PropSimPython
# Project Caelus, Aphlex 1C Engine
# Jason Chen, 10 February, 2021


import numpy as np
from scipy.integrate import solve_ivp
from typing import Tuple
from classes import Struct
from state_flow import init_liquid_state, LiquidStateVector, lsv_from_column_vec


def integration(inputs: Struct) -> Tuple[float, Struct]:
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
    record_config = { 
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
        "of_ratio_i"      : None,
        "gamma_ex"        : None,
        "m_dot_ox"        : None,
        "m_dot_fuel"      : None,
        "p_crit"          : None,
        "m_dot_ox_crit"   : None,
        "M_e"             : None,
        "p_exit"          : None,
        "p_shock"         : None
    } 
    record = Struct(record_config) # Convert into Struct class
    state_0, x0 = init_liquid_state(inputs)
    # Configure the integrator 
    solve_ivp(f, method='BDF')
    

def liquid_model(time: float, x: np.ndarray, inputs: Struct):
    """ Model engine physics for a liquid motor, provide state vector derivative. """
    
    # Constants
    R_u = inputs.constants.R_u # Universal gas constant [J/mol*K]
    M_n2o = 0.044013 # Molecular mass of nitrous oxide [kg/mol]
    R_n2o = R_u/M_n2o #  Specific gas constant of nitrous oxide [J/kg*K]
    a_n2o = 0.38828/M_n2o^2 # van der Waal's constant a for N2O [[Pa*(kg/m^3)^-2]]
    b_n2o = 44.15/M_n2o*10^-6 # van der Waal's constant a for N2O [m^3/kg]

    state = lsv_from_column_vec(x, inputs)
    state_dot = LiquidStateVector()

    # Calculate injector mass flow rate
    output = N2OTankMDot(inputs, state, time)
    m_dot_lox, m_dot_gox, m_dot_oxtank_press, T_dot_drain_ox, p_crit, m_dot_ox_crit = output


def n2o_tank_mdot(inputs: Struct, state: LiquidStateVector, time: float) -> tuple:
    """
    Calculates flow rate out of N2O Tank. Interpolates between liquid and gas flow 
    based on amount of liquid oxidizer remaining. If no liquid oxidizer remains, it
    uses gas flow. If plenty of liquid oxidizer (i.e. mass greater than set tolerance),
    liquid flow is used. Linear interpolation is used (in order to avoid hysteresis)
    in between with a small, steady liquid source in tank. 

    Inputs:
        - inputs: object representing motor inputs
        - state: object representing state of system
        - time: time
    Outputs:
        - m_dot_lox: mass flow rate out of tank of liquid oxidizer
        - m_dot_gox: mass flow rate out of tank of gaseous oxidizer
        - m_dot_oxtank_press: mass flow rate out of tank of pressurant gas
        - T_dot_drain: oxidizer tank temperature change from draining
        - p_crit: critical pressure below which draining flow is choked
        - m_dot_ox_crit: critical mass flow rate below which draining flow is choked
    """
    # Constants
    dm_lox_tol = 1e-2 # The tolerance as mentioned above for "plenty" of liquid [kg]
    M_n2o = 0.044013 # Molecular mass of nitrous oxide [kg/mol]
    a_n2o = 0.38828/M_n2o^2 # van der Waal's constant a for N2O [[Pa*(kg/m^3)^-2]]

    # TODO: Below
    """
    # If tank pressure is greater than combustion chamber pressure
    if state.p_oxmanifold > state.p_cc and inputs.ox.Cd_injector*inputs.ox.injector_area*inputs.Throttle(time)>0:
        # Find flow rates for gas, flow rates for liquid, interpolate within tolerance for smooth transition
        m_dot_lox = zeros(2,1)
        m_dot_gox = zeros(2,1)
        m_dot_oxtank_press = zeros(2,1)
        p_crit = zeros(2,1)
        m_dot_ox_crit = zeros(2,1)
        Q = zeros(2,1)
        
        ## Liquid flow
        m_dot_ox, m_dot_ox_crit(1), p_crit(1) = LN2OMDot(inputs, \
            inputs.ox.injector_area*inputs.Throttle(time), state.p_oxmanifold, state.T_oxtank, \
            state.p_cc)

        m_dot_lox(1) = m_dot_ox
        m_dot_gox(1) = 0
        m_dot_oxtank_press(1) = 0
        Q(1) = m_dot_ox/state.N2O_properties.rho_l
        
        ## Gas Flow
        d_inj = sqrt(4/pi*inputs.ox.injector_area*inputs.Throttle(time))
        [ ~, ~, ~, ~, m_dot_ox, ~ ] = NozzleCalc( d_inj, d_inj, state.T_oxtank, \
            state.p_oxmanifold, state.gamma_ox_ullage, M_n2o, state.p_cc)
        
        p_crit(2) = 0
        m_dot_ox_crit(2) = 0
        m_dot_lox(2) = 0
        m_dot_gox(2) = m_dot_ox*\
            (state.m_gox)/(state.m_gox + state.m_oxtank_press)
        m_dot_oxtank_press(2) = m_dot_ox*\
            (state.m_oxtank_press)/(state.m_gox + state.m_oxtank_press)
        Q(2) = m_dot_ox*state.V_ox_ullage/(state.m_gox + state.m_oxtank_press)
        
        ## Total Flow Rate
        if state.m_lox > dm_lox_tol:
            frac_lox = 1
        else:
            frac_lox = max(0,state.m_lox/dm_lox_tol)
        end
        
        m_dot_lox = frac_lox*m_dot_lox(1) + (1-frac_lox)*m_dot_lox(2)
        m_dot_gox = frac_lox*m_dot_gox(1) + (1-frac_lox)*m_dot_gox(2)
        m_dot_oxtank_press = frac_lox*m_dot_oxtank_press(1) + (1-frac_lox)*m_dot_oxtank_press(2)
        p_crit = frac_lox*p_crit(1) + (1-frac_lox)*p_crit(2)
        m_dot_ox_crit = frac_lox*m_dot_ox_crit(1) + (1-frac_lox)*m_dot_ox_crit(2)
        Q = frac_lox*Q(1) + (1-frac_lox)*Q(2)
    else: # bruh ur bad if ur chamber pressure is greater than tank press
        m_dot_lox = 0
        m_dot_gox = 0
        m_dot_oxtank_press = 0
        p_crit = state.p_oxtank
        m_dot_ox_crit = 0
        Q = 0
    """

# def terminal_func(t: float, x: np.ndarray, inputs: Struct):
#     # Sets the termination of integration. Integration will stop when values go to zero. 

#     state = LiquidStateVector.from_column_vector(x, inputs)
#     is_terminal = [1, 1]
#     direction = [0, 0]

#     value = [state.m_lox + state.m_gox, state.m_fuel] # Mass of propellants
#     return value, is_terminal, direction


def find_G():
    raise NotImplementedError


def find_mass_gradient():
    raise NotImplementedError
