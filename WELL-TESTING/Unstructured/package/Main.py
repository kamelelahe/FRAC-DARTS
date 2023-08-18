import pandas as pd
import numpy as np
import os
from multiprocessing import freeze_support

from .model.ModelLineSource  import LineSourceModel
from .model.ModelFracture import FractureModel
from .mesh.FractureProcessing.mesh_raw_fractures import mesh_raw_fractures

from .visualization.Visualize import Visualize
from .equations.LineSourceEq import Pwf_line_source_solution
from .equations.DimensionLess import DimensionlessTime, DimensionlessPressure

from darts.engines import redirect_darts_output
import pdb

Simulation_name = '03'
#constant input parameters
inputs = {
    'h': 1,  # m
    're': 2000,  # m
    'Ct': 1e-8,  # 1/Pa
    'mu': 1e-4,  # Pa.s
    'k':  1 ,  # mD
    'phi':0.3,
    'Pi': 200,  # bar
    'rw': 0.15,  # m
    'Q': 50,  # m3/day
    'skin' : 0,
    'Cr': 1e-10,
    'prod_well_coords': [[125, 250, 0]]
}

simulationParams = {
    'first_ts': 1e-4,
    'mult_ts': 1.2,
    'max_ts': 0.05,
    'tolerance_newton': 1e-6
}

# pdb.set_trace()
# Create a dictionary to store the results for each grid size
results = {}
size_report_step = 0.05  # Size of the reporting step (when output is writen to .vtk format)
num_report_steps =10 # Number of reporting steps (see above)
start_time = 0  # Starting time of the simulation
end_time = size_report_step * num_report_steps  # End time of the simulation

bound_cond= 'wells_in_frac'
problem_type= 'fracture' #or 'lineSource'
char_len_list = [32,64,128]
paramList=[1]
paramName='h'

# Create the ExcelWriter object outside of the loop
def main():
    for char_len in char_len_list:
        for val in paramList:
            inputs[paramName] = val
            if problem_type=='lineSource':
                # Define output directory
                outputDir ='LineSource/' + str(Simulation_name)
                # Create directories if they don't exist
                os.makedirs(outputDir, exist_ok=True)

                meshProperties = {
                    'charLength': char_len,
                    'bufferSize': 0,
                }
                m = LineSourceModel(inputs=inputs, simulationParams=simulationParams, meshProperties=meshProperties,
                          bound_cond=bound_cond, outputDir=outputDir)
            else:
                # File names and directories:
                DFN_name = '3Conjugate'
                output_dir_mesh='DFNs/' + str(DFN_name) + '/Meshes/'+'ÇharLength'+str(char_len)
                os.makedirs(output_dir_mesh, exist_ok=True)

                meshGeneration=False #this means that we want to build the mesh in place
                if meshGeneration:
                    filename_base = DFN_name
                    frac_data_raw = np.genfromtxt(os.path.join('DFNs/' + str(DFN_name) + '/Text', 'DFN001.txt'))
                    # Input parameters
                    decimals = 7  # in order to remove duplicates we need to have fixed number of decimals
                    mesh_clean = True  # need gmsh installed and callable from command line in order to mesh!!!
                    mesh_raw = True  # need gmsh installed and callable from command line in order to mesh!!!
                    num_partition_x = 4  # number of partitions for parallel implementation of intersection finding algorithm
                    num_partition_y = 4  # " ... "
                    mesh_raw_fractures(frac_data_raw, char_len, output_dir=output_dir_mesh, filename_base=filename_base,
                                       height_res=1, apertures_raw=None,
                                       box_data=None, margin=25, mesh_raw=mesh_raw, decimals=decimals,
                                       tolerance_zero=1e-10, tolerance_intersect=1e-10,
                                       calc_intersections_before=True, num_partition_x=num_partition_x,
                                       num_partition_y=num_partition_y,
                                       partition_fractures_in_segms=True, matrix_perm=1, char_len_mult=1,
                                       char_len_boundary=None, main_algo_iters=1)

                mesh_name = DFN_name + '_raw_lc_' + str(char_len) + '.msh'
                # Define output directory
                mesh_path = output_dir_mesh + '/' + mesh_name
                m = FractureModel(inputs=inputs, simulationParams=simulationParams,
                          bound_cond=bound_cond, mesh_path=mesh_path)



                # Define output directory
                outputDir ='DFNs/' + str(DFN_name)+'/results/' + str(Simulation_name)
                # Create directories if they don't exist
                os.makedirs(outputDir, exist_ok=True)

                # redirect output
                redirect_darts_output(file_name=str(Simulation_name) + '.log')

            print('============================================================================================')
            print('++++++++++++++++++++++++++' + str(char_len) + '+++++++++++++++++++++++++++++++++++++')
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


            with pd.ExcelWriter(f'{outputDir}/{Simulation_name}.xlsx', engine='openpyxl') as writer:

                for ith_step in range(num_report_steps):

                    # Create the model with the current grid size and refinement
                    print('================= Unstructured === Step :'+ str(ith_step+1) + '======================================')
                    m.run_python(size_report_step)

                    # Store the results for the current grid size and refinement
                    SimName= Simulation_name+'char_len'+str(char_len)
                    data = pd.DataFrame.from_dict(m.physics.engine.time_data)
                    results[(char_len,val)] = data

                    # Use the writer to write to a different sheet
                    sheet_name='char_len'+str(char_len)
                    data.to_excel(writer, sheet_name=sheet_name)

                    # Prepare property_array and cell_property
                    pressure_field = m.physics.engine.X[:-1:2]
                    saturation_field = m.physics.engine.X[1::2]
                    property_array = np.empty((tot_unknws, tot_properties))
                    property_array[:, 0] = pressure_field
                    property_array[:, 1] = saturation_field

                    # Write to VTK
                    outputDir_VTK = outputDir+'/vtk_data'
                    if not os.path.exists(outputDir_VTK):
                        os.makedirs(outputDir_VTK)
                    fileName = 'solution' + str(char_len)
                    m.reservoir.unstr_discr.write_to_vtk(outputDir_VTK, property_array, m.physics.vars, ith_step,fileName)

                    # Print timers and statistics for the run
                    m.print_timers()
                    m.print_stat()

                    #generating visualized outputs
                    v = Visualize(results=results,paramList=paramList,paramName=str(paramName), labels=['CharacteristicLength'], output_dir=outputDir, inputs=inputs)
                    frameName = 'solution'+str(char_len) + str(ith_step)
                    #v.showVtk(fileName=frameName)
                    v.imageVtk(fileName=frameName,OuutputName=frameName+'Grid',output_dir=outputDir, showImage=True, gridding=True)
                    v.imageVtk(fileName=frameName,OuutputName=frameName+'3D',output_dir=outputDir, showImage=True, gridding=False)
                    v.image2DVtk(fileName=frameName,OuutputName=frameName,output_dir=outputDir, showImage=True, gridding=False)
                v.makeGif(input_name='solution'+str(char_len), num_report_steps=num_report_steps, output_file=SimName + '.gif', duration=2)

    ## visualization ##
    v = Visualize(results=results,paramList=paramList,paramName=str(paramName), labels=['CharacteristicLength',str(paramName)], output_dir=outputDir, inputs=inputs)

    if problem_type=='lineSource':

         ## Calculating analytical solution ##
        deltaT_analytical= simulationParams['first_ts']
        t_values = np.arange(simulationParams['first_ts'], size_report_step*num_report_steps, deltaT_analytical)

        # Initialize an empty dictionary to store Pwf values
        Pwf_values_dic = {}
        # Loop over t_values and calculate Pwf for each t
        for val in paramList:
            Pwf_values_list = []
            for t in t_values:
                inputs[paramName]= val
                Pwf = Pwf_line_source_solution(inputs,t)
                Pwf_values_list.append(Pwf)
            Pwf_values_dic[val] = Pwf_values_list

        ## Conversion to DimensionLess Parameters ##
        tD_analytical = {}
        pD_analytical_dic = {}
        for val in paramList:
            inputs[paramName]= val
            tD_analytical[val] = DimensionlessTime(inputs, t_values)
            pD_analytical_dic[val] = DimensionlessPressure(inputs, Pwf_values_dic[val])

        # plotting with analytical solution
        v.plot_pD_tD(fileName='pD_tD_ana', pD_analytical=pD_analytical_dic, tD_analytical=tD_analytical, plotAnalytical=True, savePlot=True, showPlot=False)
        v.plot_pD_tD_log(fileName= 'pD_tD_log_ana',pD_analytical=pD_analytical_dic, tD_analytical=tD_analytical, plotAnalytical=True,savePlot=True, showPlot=False)
        v.plot_pressureVsTime(fileName='P_T_ana' , p_analytical=Pwf_values_dic,t_analytical=t_values, plotAnalytical=True, savePlot=True, showPlot=False)
        v.plot_pressureVsTime_semiLog(fileName='p_t_log_ana' ,p_analytical=Pwf_values_dic,t_analytical=t_values, plotAnalytical=True, savePlot=True, showPlot=False)

    else:

        # plotting whiout analytical solution
        v.plot_pD_tD(fileName='pD_tD', plotAnalytical=False, savePlot=True, showPlot=False)
        v.plot_pD_tD_log(fileName='pD_tD_log', plotAnalytical=False, savePlot=True, showPlot=False)
        v.plot_pressureVsTime(fileName='P_T', plotAnalytical=False, savePlot=True, showPlot=False)
        v.plot_pressureVsTime_semiLog(fileName='P_T_log', plotAnalytical=False, savePlot=True, showPlot=False)
if __name__ == "__main__":
    freeze_support()
    main()
