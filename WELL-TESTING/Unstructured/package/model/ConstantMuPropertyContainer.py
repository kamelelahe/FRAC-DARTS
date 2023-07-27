from darts.physics.properties.basic import ConstFunc
from darts.engines import value_vector
from darts.physics.properties.iapws.iapws_property import *
from darts.physics.properties.iapws.custom_rock_property import *

class MyPropertyContainer:
    '''
    Class responsible for collecting all needed properties in geothermal simulation
    '''
    def __init__(self, inputs):
        """
        here we modified property container to accept some values from input parameters
        """
        temp = 350
        self.rock = [value_vector([inputs['Pi'], inputs['Ct'], temp])]

        # properties implemented in python (the IAPWS package)
        self.temperature = ConstFunc(temp)      # Create temperature object
        self.water_enthalpy = iapws_water_enthalpy_evaluator()  # Create water_enthalpy object
        self.steam_enthalpy = iapws_steam_enthalpy_evaluator()  # Create steam_enthalpy object
        self.water_saturation = ConstFunc(1)  # Create water_saturation object  with constant value
        self.steam_saturation = ConstFunc(0)  # Create steam_saturation object with constant value
        self.water_relperm = ConstFunc(1)   # Create water_relperm object with constant value
        self.steam_relperm = ConstFunc(0)   # Create steam_relperm object with constant value
        self.water_density = iapws_water_density_evaluator()  # Create water_density object
        self.steam_density = iapws_steam_density_evaluator()  # Create steam_density object
        self.water_viscosity = ConstFunc(inputs['mu']*1e3)    # Create water_viscosity object with constant value
        self.steam_viscosity = iapws_steam_viscosity_evaluator()          # Create steam_viscosity object
        self.rock_compaction = custom_rock_compaction_evaluator(self.rock)  # Create rock_compaction object
        self.rock_energy = custom_rock_energy_evaluator(self.rock)          # Create rock_energy object


