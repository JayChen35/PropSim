# Prop Sim Classes
# Project Caelus 2/4/2021


class Inputs():
    """
    Stores inputs for propsim
    """

    def __init__(self):
        self.ox = None
        self.ox_pressurant = None
        
        self.fuel = None
        self.fuel_pressurant = None
        
        self.dt_valve_open = None
        self.mass_dry_rocket = None
        
        self.length_cc = None
        self.d_cc = None
        
        self.nozzle_efficiency = None
        self.nozzle_correction_factor = None
        self.c_star_efficiency = None
        
        self.d_throat = None
        self.exp_ratio = None
        
        self.T_amb = None
        self.p_amb = None
        
        self.comb_data = None
        
        
        #ox: injector_area, Cd_injector , V_l, V_tank, tank_id, h_offset_tank, d_flowline, T_tank
        # fuel: injector_area, Cd_injector, h_offset_tank, d_flowline, rho
        # dt_valve_open
        # mass_dry_rocket
        # d_cc
        # length_cc
        
        # ox_pressurant: gas_properties, set_pressure, storage_initial_pressure, tank_volume, flow_CdA
        # fuel_pressurant: gas_properties, set_pressure, storage_initial_pressure, tank_volume, flow_CdA
        

class Gas():
    """
    Stores common information about gasses
    """ 
    
    def __init__(self):
        self.c_v = None
        self.molecular_mass = None
        
class Pressurant():
    """
    PRESSURANT Represents a propellant pressurant mechanism
    Represents an instance of a propellant under pressurization by a
    gas.
    """
    
    def __init__(propellant_type = None):
        # Whether pressurant is to be simulated
        self.gas_properties = None
        self.set_pressure = None
        self.storage_initial_pressure = None
        self.tank_volume = None
        self.flow_CdA = None
        if propellant_type == 'fuel' or propellant_type == 'oxidizer':
            self.propellant_type = propellant_type
        else:
            raise Exception("Invalid Propellant Type")