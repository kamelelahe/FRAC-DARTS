import pandas as pd
import numpy as np
import os
from .model.Model  import Model
from .visualization.Visualize import Visualize
from .equations.LineSourceEq import Pwf_line_source_solution
from .equations.DimensionLess import DimensionlessTime, DimensionlessPressure
from darts.engines import redirect_darts_output
Simulation_name= '04Unstructured'
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
    'Q': 20,  # m3/day
    'skin':0
}

simulationParams = {
    'first_ts': 1e-3,
    'mult_ts': 4,
    'max_ts': 0.001,
    'tolerance_newton': 1e-4
}

# List of grid sizes
grid_sizes = [50]
# List of grid refinements
grid_refinements = [4]
#list of grid refinement in z direction
Dz = [4]
# Create a dictionary to store the results for each grid size
results = {}
simulation_time = 0.1
MuList=[1e-3]
#MuList= [1e-3, 7e-4, 3774e-4, 1e-4]
# Create the ExcelWriter object outside of the loop

with pd.ExcelWriter(f'{outputDir}/{Simulation_name}.xlsx', engine='openpyxl') as writer:
    # Loop over the grid sizes
    for grid_size in grid_sizes:
        # Loop over the grid refinements
        for grid_refinement in grid_refinements:
            #loop over size of grid in z durection
            for dz in Dz:
                for Mu in MuList:
                    inputs['mu']= Mu

                    # Create the model with the current grid size and refinement
                    print('============================================================================================')
                    print('========'+' Mu= '+str(Mu)+'GS='+str(grid_size)+'GR='+str(grid_refinement)+'Dz='+str(dz)+str(Mu))
                    m = Model(inputs=inputs, simulationParams=simulationParams, dx=grid_size, dy=grid_size, dz=dz,  grid_refinement=grid_refinement)
                    m.init()

                    # Run the simulation
                    m.run_python(simulation_time)

                    # Print timers and statistics for the run
                    m.print_timers()
                    m.print_stat()
                    # Store the results for the current grid size and refinement
                    SimName='GS'+str(grid_size)+'GR'+str(grid_refinement)+'Dz'+str(dz)+str(Mu)
                    data = pd.DataFrame.from_dict(m.physics.engine.time_data)
                    results[(grid_size, grid_refinement, dz,Mu,)] = data
                    # Use the writer to write to a different sheet
                    data.to_excel(writer, sheet_name=SimName)
                    original_dir = os.getcwd()
                    os.chdir(outputDir)
                    m.export_pro_vtk(file_name=SimName)
                    os.chdir(original_dir)

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



## visualization ##
v = Visualize(results=results,MuList=MuList, labels=['GS', 'GR', 'dz','Mu'], output_dir=outputDir, inputs=inputs)

"""
# plotting whiout analytical solution
v.plot_pD_tD(fileName='pD_tD_ana', plotAnalytical=False, savePlot=True, showPlot=False)
v.plot_pD_tD_log(fileName= 'pD_tD_log_ana', plotAnalytical=False,savePlot=True, showPlot=False)
v.plot_pressureVsTime(fileName='P_T_ana' , plotAnalytical=False, savePlot=True, showPlot=False)
v.plot_pressureVsTime_semiLog(fileName='p_t_log_ana' , plotAnalytical=False, savePlot=True, showPlot=False)
"""

# plotting with analytical solution
v.plot_pD_tD(fileName='pD_tD_ana', pD_analytical=pD_analytical_dic, tD_analytical=tD_analytical, plotAnalytical=True, savePlot=True, showPlot=False)
v.plot_pD_tD_log(fileName= 'pD_tD_log_ana',pD_analytical=pD_analytical_dic, tD_analytical=tD_analytical, plotAnalytical=True,savePlot=True, showPlot=False)
v.plot_pressureVsTime(fileName='P_T_ana' , p_analytical=Pwf_values_dic,t_analytical=t_values, plotAnalytical=True, savePlot=True, showPlot=False)
v.plot_pressureVsTime_semiLog(fileName='p_t_log_ana' ,p_analytical=Pwf_values_dic,t_analytical=t_values, plotAnalytical=True, savePlot=True, showPlot=False)

for grid_size in grid_sizes:
    for grid_refinement in grid_refinements:
        for dz in Dz:
            for Mu in MuList:
                frameName='GS'+str(grid_size)+'GR'+str(grid_refinement)+'Dz'+str(dz)+str(Mu)+'_ts0'
                #v.showVtk(fileName=frameName)
                v.imageVtk(fileName=frameName, output_dir=outputDir, showImage=True )


