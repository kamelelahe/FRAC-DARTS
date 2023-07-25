import numpy as np
import math
from .ConstantMuPropertyContainer import MyPropertyContainer
from darts.models.reservoirs.struct_reservoir import StructReservoir
from darts.physics.geothermal.physics import  Geothermal
from darts.models.darts_model import DartsModel

class Model(DartsModel):

    def __init__(self, inputs, simulationParams, dx, dy, dz, grid_refinement=1, n_points=128):
        """in init part, we initialize reservoir parameters. static properties, grid properties, location of wells, and petropysical proprties are add
        this stage"""
        # Call base class constructor
        super().__init__()

        # Start the timer for initialization
        self.timer.node["initialization"].start()

        # Extract inputs
        self.flow_rate = inputs['Q']
        viscosity = inputs['mu']
        permeability = inputs['k']
        reservoir_height = inputs['h']
        reservoir_radius = inputs['re']
        self.well_bore_radius = inputs['rw']
        self.initial_pressure = inputs['Pi']
        compressibility = inputs['Ct']
        self.skin = inputs['skin']

        time = 0  # Time in hours

        # calculate equivalent square side length based on radus of reservoir
        reservoir_width = math.sqrt(math.pi) * reservoir_radius
        """ Adjust the gridding parameters"""
        self.grid_refinement = grid_refinement
        dx_large = dx
        dy_large = dy
        self.dz = dz
        # large blocks
        (nx_large, ny_large, self.nz) = (round(reservoir_width / dx_large) - 3, round(reservoir_width / dy_large) - 3,
                                         round(
                                             reservoir_height / self.dz))  # Number of grid blocks in x, y, z directions
        (nx_small, ny_small) = (
        3 * self.grid_refinement, 3 * self.grid_refinement)  # Number of grid blocks in x, y, z directions
        self.nx = nx_large + nx_small
        self.ny = ny_large + ny_small
        self.nx_large = nx_large
        self.nx_small = nx_small
        # grid refinement
        self.nb_Total = self.nx * self.ny * self.nz

        """ set up location of well"""
        # well_x, well_y = (round(reservoir_width/2),round(reservoir_width/2))
        # Define well's locations
        self.jw = [round(self.nx / 2)]
        self.iw = [round(self.ny / 2)]

        """ Set-up reservoir properties"""
        # Define parameters for the reservoir

        perm = np.ones(self.nb_Total) * permeability  # Create perm array filled with 100
        poro = np.ones(self.nb_Total) * 0.3  # Create porosity array filled with 0.3

        self.dx = np.ones(self.nb_Total) * dx_large
        coordinate = 0 - self.nx
        for k in range(self.nz):
            for j in range(self.ny):
                coordinate = coordinate + self.nx
                for i in range(self.iw[0] - round(1.5 * self.grid_refinement),
                               self.iw[0] + round(1.5 * self.grid_refinement)):
                    self.dx[coordinate + i] = dx_large / self.grid_refinement

        self.dy = np.ones(self.nb_Total) * dy_large
        for k in range(self.nz):
            for j in range(self.jw[0] - round(1.5 * self.grid_refinement),
                           self.jw[0] + round(1.5 * self.grid_refinement)):
                self.dy[
                j * self.nx + k * self.nx * self.ny: j * self.nx + self.nx + k * self.nx * self.ny] = dy_large / self.grid_refinement

        # Create structured reservoir with defined properties
        self.reservoir = StructReservoir(self.timer, nx=self.nx, ny=self.ny, nz=self.nz, dx=self.dx, dy=self.dy,
                                         dz=self.dz, permx=perm,
                                         permy=perm, permz=perm, poro=poro, depth=2000)

        # Set open boundaries
        self.reservoir.set_boundary_volume(xz_minus=1e8, xz_plus=1e8, yz_minus=1e8, yz_plus=1e8, xy_minus=1e6)

        """ Specifically for geothermal reservoir"""
        # Set rock heat capacity and rock thermal conduction
        hcap = np.array(self.reservoir.mesh.heat_capacity, copy=False)
        rcond = np.array(self.reservoir.mesh.rock_cond, copy=False)
        hcap.fill(2200)

        #self.property_container = PropertyContainer()
        self.property_container= MyPropertyContainer(inputs=inputs)
        # create pre-defined physics for geothermal
        self.physics = Geothermal(self.timer, property_container=self.property_container,
                                  n_points=n_points, min_p=1, max_p=351, min_e=1000, max_e=10000, cache=False)
        self.physics.init_physics()

        """ adding wells and constraints"""

        # Add production well
        self.reservoir.add_well("PRD")
        n_perf = self.nz
        # Add perforations to the payzone
        for n in range(1, n_perf):
            self.reservoir.add_perforation(well=self.reservoir.wells[-1], i=self.iw[0], j=self.jw[0], k=n + 1,
                                           well_radius=self.well_bore_radius)

        """ Linear solver"""
        # Set timestep parameters
        self.params.first_ts = simulationParams['first_ts']
        self.params.mult_ts = simulationParams['mult_ts']
        self.params.max_ts = simulationParams['max_ts']
        self.params.tolerance_newton = simulationParams['tolerance_newton']

        # Stop the timer for initialization
        self.timer.node["initialization"].stop()

    def set_initial_conditions(self):
        """initial pressure and temprature"""
        # Set uniform initial conditions for pressure and temperature
        self.physics.set_uniform_initial_conditions(self.reservoir.mesh, uniform_pressure=self.initial_pressure,
                                                    uniform_temperature=350)

    def set_boundary_conditions(self):
        """ boundry condition"""
        # Set boundary conditions for wells
        # For the producer, we control the rate
        for i, w in enumerate(self.reservoir.wells):
            w.control = self.physics.new_rate_water_prod(self.flow_rate)

    def export_pro_vtk(self, file_name='Results', output_dir=''):
        # connect to simulation array
        X = np.array(self.physics.engine.X, copy=False)
        nb = self.reservoir.mesh.n_res_blocks
        # compute temperature using pressure and enthalpy (in different units)
        #temp = _Backward1_T_Ph_vec(X[0:2 * nb:2] / 10, X[1:2 * nb:2] / 18.015)
        # define additional arrays to the output
        local_cell_data = {'Perm': self.reservoir.global_data['permx']}
        # use export_vtk defined in the base class (DartsModel)
        self.export_vtk(file_name=file_name, local_cell_data=local_cell_data)
