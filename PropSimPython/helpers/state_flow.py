# Helper Methods for PropSimPython
# Project Caelus, Aphlex 1C Engine
# 02 February, 2021

import numpy as np


class LiquidStateVec():
    def __init__(self, inputs: dict):
        super().__init__()
        self.T_fuel_tank_press = None # Fuel tank temperature
        self.m_fuel_tank_press = None # Mass of pressurant in oxidizer tank
        self.m_fuel_press_tank = None # Mass of pressurant in oxidizer pressurant tank
        self.m_fuel = None # Fuel mass, kg
        self.inputs = inputs
    
    def p_fuel_tank(self):
        # Calculate fuel tank pressure
        p_fueltank = self.m_fuel_tank_press*self.inputs.fuel_pressurant.gas_properties.R_specific* \
            self.T_fueltank_press/self.V_fuel_ullage
        return p_fueltank
    
    # def p_fuelpresstank(self):
    #     # Calculate pressure in fuel pressurant tank
    #     m_press_initial = self.inputs.fuel_pressurant.storage_initial_pressure* \
    #         self.inputs.fuel_pressurant.tank_volume/ \
    #             (self.inputs.fuel_pressurant.gas_properties.R_specific.*self.inputs.T_amb)
    #     p_fuelpresstank = self.inputs.fuel_pressurant.storage_initial_pressure* \
    #         (self.m_fuelpresstank/m_press_initial)^self.inputs.fuel_pressurant.gas_properties.gamma
    #     return p_fuelpresstank
    
    # def T_fuelpresstank = T_fuelpresstank(self):
    #     # Calculate temperature in oxidizer pressurant tank
    #     T_fuelpresstank = self.inputs.T_amb*
    #         (self.p_fuelpresstank/self.inputs.fuel_pressurant.storage_initial_pressure)^
    #         ((self.inputs.fuel_pressurant.gas_properties.gamma-1)/
    #         self.inputs.fuel_pressurant.gas_properties.gamma);
    
    def V_fuel(self):
        # Calculate volume of fuel
        V_fuel = self.m_fuel / self.inputs["fuel"]["rho"]
        return V_fuel
    
    def V_fuel_ullage(self):
        # Calculate ullage volume in fuel tank
        V_fuel_ullage = self.inputs["fuel"]["tank_V"] - self.V_fuel()
        return V_fuel_ullage
    
    def V_cc(self):
        # Calculate combustion chamber volume
        V_cc = np.pi/4*(self.inputs["d_cc"])**2*self.inputs["length_cc"]
        return V_cc
    
    def column_vec(self):
        # Create output column vector
        column_vector = [self.m_lox;
            self.m_gox;
            self.m_oxtank_press;
            self.m_oxpresstank;
            self.m_fuel_tank_press;
            self.m_fuelpresstank;
            self.T_oxtank;
            self.T_fueltank_press;
            self.m_fuel;
            self.m_cc;
            self.M_cc;
            self.gamma_cc;
            self.T_cc;];
        return column_vector


def init_liquid_state(inputs: dict) -> Tuple[LiquidStateVec, np.ndarray]:
    # Initializes the state vector for a liquid system.
    # Uses the inputs to create an initial state vector for a liquid
    state_0 = LiquidStateVector(inputs);
    state_0 = InitializeOxtank(state_0, inputs);
    state_0 = InitializeFueltank(state_0, inputs);
    state_0 = InitializeCombustionChamber(state_0, inputs);
    x0 = state_0.ColumnVector;
    return [state_0, x0]
