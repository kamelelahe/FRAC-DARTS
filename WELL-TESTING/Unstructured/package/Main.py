import pandas as pd
import numpy as np
import os
import meshio
from .model.Model  import Model
from .visualization.Visualize import Visualize
from .equations.LineSourceEq import Pwf_line_source_solution
from .equations.DimensionLess import DimensionlessTime, DimensionlessPressure
from darts.engines import redirect_darts_output

Simulation_name = '15UnStructured-NewMesh'
outputDir = 'results/'+str(Simulation_name)
if not os.path.exists(outputDir):
    os.makedirs(outputDir)
redirect_darts_output(file_name=str(Simulation_name)+'.log')

#constant input parameters
inputs = {
    'h': 40,  # m
    're': 2000,  # m
    'Ct': 1e-8,  # 1/Pa
    'mu': 0.3774e-3,  # Pa.s (viscousity of water at 350K)
    'k':  2533.22 ,  # mD
    'Pi': 300,  # bar
    'rw': 0.15,  # m
    'Q': 400,  # m3/day
    'skin':0
}

simulationParams = {
    'first_ts': 0.5,
    'mult_ts': 4,
    'max_ts': 0.5,
    'tolerance_newton': 1e-4
}

# Create a dictionary to store the results for each grid size
results = {}
size_report_step = 1  # Size of the reporting step (when output is writen to .vtk format)
num_report_steps = 10 # Number of reporting steps (see above)
start_time = 0  # Starting time of the simulation
end_time = size_report_step * num_report_steps  # End time of the simulation
bound_cond= 'const_pres_rate'
#mesh_file='C:\meshes\mesh_60_real_6.msh'
mesh_file='C:\meshes\Ilshat.msh'
###initialization ####
m = Model(inputs=inputs, simulationParams=simulationParams,
          mesh_file=mesh_file, bound_cond=bound_cond)
m.init()

# Prepare property_array and cell_property
tot_unknws = m.reservoir.unstr_discr.fracture_cell_count + m.reservoir.unstr_discr.matrix_cell_count + len(
    m.reservoir.wells) * 2
tot_properties = 2
pressure_field = m.physics.engine.X[:-1:2]
saturation_field = m.physics.engine.X[1::2]
property_array = np.empty((tot_unknws, tot_properties))
property_array[:, 0] = pressure_field
property_array[:, 1] = saturation_field

# Write to VTK
#m.reservoir.unstr_discr.write_to_vtk(outputDir, property_array, m.physics.vars, 0)

# Create the ExcelWriter object outside of the loop
with pd.ExcelWriter(f'{outputDir}/{Simulation_name}.xlsx', engine='openpyxl') as writer:
    for ith_step in range(num_report_steps):
        # Create the model with the current grid size and refinement
        print('============================================================================================')
        print('================= Unstructured === Step:'+ str(ith_step) + '======================================')
        m.run_python(size_report_step)
        # Store the results for the current grid size and refinement
        SimName=Simulation_name
        data = pd.DataFrame.from_dict(m.physics.engine.time_data)
        results[SimName] = data
        # Use the writer to write to a different sheet
        data.to_excel(writer, sheet_name=SimName)
        # Prepare property_array and cell_property
        pressure_field = m.physics.engine.X[:-1:2]
        saturation_field = m.physics.engine.X[1::2]
        property_array = np.empty((tot_unknws, tot_properties))
        property_array[:, 0] = pressure_field
        property_array[:, 1] = saturation_field
        # Write to VTK
        m.reservoir.unstr_discr.write_to_vtk(outputDir, property_array, m.physics.vars, ith_step)
        # Print timers and statistics for the run
        m.print_timers()
        m.print_stat()

## visualization ##
v = Visualize(results=results, labels=['SimName'], output_dir=outputDir, inputs=inputs)


# plotting whiout analytical solution
v.plot_pD_tD(fileName='pD_tD', plotAnalytical=False, savePlot=True, showPlot=False)
v.plot_pD_tD_log(fileName= 'pD_tD_log', plotAnalytical=False,savePlot=True, showPlot=False)
v.plot_pressureVsTime(fileName='P_T', plotAnalytical=False, savePlot=True, showPlot=False)
v.plot_pressureVsTime_semiLog(fileName='p_t_log', plotAnalytical=False, savePlot=True, showPlot=False)

for ith_step in range(num_report_steps):
    frameName='solution'+str(ith_step)
    #v.showVtk(fileName=frameName)
    v.imageVtk(fileName=frameName, pressure_range=(200 ,300))

#v.makeGif( num_report_steps=num_report_steps, output_file=Simulation_name+'.gif', duration=2)

"""
###Calculation for numerical solution###
## Calculating analytical solution ##
deltaT_analytical= simulationParams['max_ts']
t_values = np.arange(simulationParams['first_ts']+ deltaT_analytical, simulation_time, deltaT_analytical)

# Initialize an empty dictionary to store Pwf values
Pwf_values_dic = {}
# Loop over t_values and calculate Pwf for each t
for Mu in MuList:
    Pwf_values_list = []
    for t in t_values:
        inputs['mu']=Mu
        Pwf = Pwf_line_source_solution(inputs,t)
        Pwf_values_list.append(Pwf)
    Pwf_values_dic[Mu] = Pwf_values_list

## Conversion to DimensionLess Parameters ##
tD_analytical = {}
pD_analytical_dic = {}
for Mu in MuList:
    inputs['mu']=Mu
    tD_analytical[Mu] = DimensionlessTime(inputs, t_values)
    pD_analytical_dic[Mu] = DimensionlessPressure(inputs, Pwf_values_dic[Mu])
    
# plotting with analytical solution
v.plot_pD_tD(fileName='pD_tD_ana', pD_analytical=pD_analytical_dic, tD_analytical=tD_analytical, plotAnalytical=True, savePlot=True, showPlot=False)
v.plot_pD_tD_log(fileName= 'pD_tD_log_ana',pD_analytical=pD_analytical_dic, tD_analytical=tD_analytical, plotAnalytical=True,savePlot=True, showPlot=False)
v.plot_pressureVsTime(fileName='P_T_ana' , p_analytical=Pwf_values_dic,t_analytical=t_values, plotAnalytical=True, savePlot=True, showPlot=False)
v.plot_pressureVsTime_semiLog(fileName='p_t_log_ana' ,p_analytical=Pwf_values_dic,t_analytical=t_values, plotAnalytical=True, savePlot=True, showPlot=False)
"""



