# Helper Methods for PropSimPython
# Project Caelus, Aphlex 1C Engine
# 02 February, 2021


import numpy as np


def n2o_properties(temp: int or float) -> dict:
    """
    Calculates an array of properties of nitrous oxide given a temperature in K
    WARNING: if temperature input is outside of -90 to 30 C, properties will
    be generated for boundary (-90 C or 30 C)
    """
    properties = dict() #creates properties (output array) as a dictionary for which data from text file can be entered
    properties["Pvap"]       = None # nitrous vapor pressure
    properties["rho_l"]      = None # liquid density of nitrous
    properties["rho_g"]      = None # gas denisty of nitrous
    properties["deltaH_vap"] = None #
    properties["cp_l"]       = None # ?
    properties["cv_l"]       = None # specific volume of liquid nitrous
    properties["cp_g"]       = None # 
    properties["cv_g"]       = None # specific volume of gaseous nitrous 
    properties["h_l"]        = None # enthalpy of liquid nitrous
    properties["h_g"]        = None # enthalpy of gaseous nitrous
    properties["s_l"]        = None # entropy of liquid nitrous
    properties["s_g"]        = None # entropy of gaseous nitrous
    properties["mu_l"]       = None # dynamic viscosity of nitrous liquid
    properties["mu_g"]       = None # dynamic viscosity of nitrous gas
    properties["e_l"]        = None # 
    properties["e_g"]        = None
        
    R_u = 8.3144621 # Universal gas constant [J/mol*K]
    M_n2o = 0.044013 # Molecular mass of nitrous oxide [kg/mol]
    R_n2o_0 = R_u/M_n2o # Specific gas constant of nitrous oxide [J/kg*K]

    # Range-check temperature
    if temp < (-90 + 273.15): # if temperature is less that bottom bound
        temp = -90 + 273.150001  # temperature is bottom bound for interpolate
    elif temp > (30 + 273.150001): # if temperature greater than top bound, 
        temp = 30 + 273.150001 # temperature equal to top bound for interpolate

    Tcrit = 309.57 #K
    Pcrit = 7251   #kPa
    rhocrit = 452 #kg/m^3
    #possibly add critical compressibility factor "Z"
    Tr = temp/Tcrit

    # Calculate vapor pressure, valid -90 to 36C
    b1 = -6.71893
    b2 = 1.3596
    b3 = -1.3779
    b4 = -4.051

    properties["Pvap"] = np.exp((1/Tr)*(b1*(1-Tr) + b2*(1-Tr)**(3/2) + b3*(1-Tr)**(5/2) + b4*(1-Tr)**5))*Pcrit
    properties["Pvap"] = properties["Pvap"]*1000

    # Calculate Density of Liquid, valid -90C to 36C
    b1 = 1.72328
    b2 = -0.83950
    b3 = 0.51060
    b4 = -0.10412

    properties["rho_l"] = np.exp(b1*(1-Tr)**(1/3) + b2*(1-Tr)**(2/3) + b3*(1-Tr) + b4*(1-Tr)**(4/3))*rhocrit

    # Calculate Density of Gas, valid -90C to 36C
    b1 = -1.00900
    b2 = -6.28792
    b3 = 7.50332
    b4 = -7.90463
    b5 = 0.629427
    Trinv = 1./Tr

    properties["rho_g"] = np.exp(b1*(Trinv-1)**(1/3) + b2*(Trinv-1)**(2/3) + b3*(Trinv-1) + b4*(Trinv-1)**(4/3) + b5*(Trinv-1)**(5/3))*rhocrit

    # Calculate dynamic viscosity of saturated liquid, valid from -90C to 30C
    b1 = 1.6089
    b2 = 2.0439
    b3 = 5.24
    b4 = 0.0293423
    theta = (Tcrit-b3)/(temp-b3)

    properties["mu_l"] = b4*np.exp(b1*(theta-1)**(1/3) + b2*(theta-1)**(4/3))

    # Calculate dynamic viscosity of saturated vapor, valid from -90C to 30C
    b1 = 3.3281
    b2 = -1.18237
    b3 = -0.055155
    Trinv = 1./Tr

    properties["mu_g"] = np.exp(b1 + b2*(Trinv-1)**(1/3) + b3*(Trinv-1)**(4/3))

    # TODO: Find Specific Enthalpy?

    NIST_data = dict() # NIST_data is an array that stores variables regarding nitrous
    
    reader = open("N2O_Properties.cgi.txt", 'r') # create the dictonary (reader) that takes information from N2O_Properties.cgi.txt
    try:
        tempList = reader.readline().split("\t")
        arr = [list() for x in range(len(tempList))]
        for x in range(0,126):
            tempList = reader.readline().split("\t") # read each line in the N2O_Properties.cgi.txt document and enter each line into the dictionary, separated by tabs
            for i, val in enumerate(tempList):
                arr[i].append(val)
        NIST_data["T"] =     arr[0]
        NIST_data["h_liq"] = arr[5]
        NIST_data["h_gas"] = arr[17]
        NIST_data["e_liq"] = arr[4]
        NIST_data["e_gas"] = arr[16]
        NIST_data["cv_l"]  = arr[7]
        NIST_data["cv_g"]  = arr[19]
        NIST_data["s_l"]   = arr[6]
        NIST_data["s_g"]   = arr[18]

    reader.close()
    
    #properties["h_l"] = np.interp(NIST_data["T"], NIST_data.h_liq, temp)*1000 # J/kg, 

    # Gas Specific Enthalpy
    properties["h_l"] = np.interp(temp, NIST_data["T"], NIST_data["h_liq"])*1000 # J/kg, 
    properties["h_g"] = np.interp(temp, NIST_data["T"], NIST_data.["h_g"])*1000 # J/kg
    properties["e_l"] = np.interp(temp, NIST_data["T"], NIST_data.["e_liq"])*1000 # J/kg
    properties["e_g"] = np.interp(temp, NIST_data["T"], NIST_data.["e_gas"])*1000 # J/kg
    properties["deltaH_vap"] = properties["h_g"]-properties["h_l"]
    properties["deltaE_vap"] = properties["e_g"]-properties["e_l"]

    #np.interp(temp, )

    # Specific Heat at Constant Volume
    properties["cv_l"] = np.interp(temp, NIST_data["T"], NIST_data["cv_l"])*1000
    properties["cv_g"] = np.interp(temp, NIST_data["T"], NIST_data["cv_g"])*1000

    # Specific Heat at Constant Pressure
    properties["cp_l"] = np.interp(temp, NIST_data["T"], NIST_data["cp_l"])*1000
    properties["cp_g"] = np.interp(temp, NIST_data["T"], NIST_data["cp_g"])*1000

    # Specific Entropy
    properties["s_l"] = np.interp(temp, NIST_data["T"], NIST_data["s_l"])*1000
    properties["s_g"] = np.interp(temp, NIST_data["T"], NIST_data["s_g"])*1000

    # Convert Properties to Standard Units
    properties["mu_l"] = properties["mu_l"]*10^-3 # mN*s/(m^2) -> N*s/m^2
    properties["mu_g"] = properties["mu_g"]*10^-6 # uN*s/(m^ 2)-> N*s/m^2

    return properties

def n2o_find_T(p_vap: int or float) -> float:
    raise NotImplementedError