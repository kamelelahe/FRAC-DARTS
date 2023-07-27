import matplotlib.pyplot as plt
import numpy as np
from ..equations.DimensionLess import DimensionlessPressure, DimensionlessTime
import matplotlib.cm as cm
import pyvista as pv
from IPython.display import Image
from PIL import Image as PILImage
import os
import imageio.v2 as imageio
import os
class Visualize:

    def __init__(self, results, inputs, labels, output_dir):
        self.output_dir = output_dir
        self.results = results
        self.inputs = inputs
        self.labels = labels
        # Create a colormap
        self.colors = cm.rainbow(np.linspace(0, 1, len(self.results)))

    def plot_pressureVsTime(self, p_analytical=[], t_analytical=[] ,plotAnalytical=True,  showPlot=True, savePlot=False, fileName='P_T'):
        # Plot the results
        fig, ax = plt.subplots(figsize=(11,8))
        # Loop over the results
        for (color, (keys, df)) in zip(self.colors, self.results.items()):
            # Plot data
            label = ', '.join(f'{label}: {value}' for label, value in zip(self.labels, keys))
            ax.plot(df['time'], df['PRD : BHP (bar)'], marker='o', linestyle='-', label=label, color=color)

        if plotAnalytical:
            ax.plot(t_analytical, p_analytical, color='r')

        ax.tick_params(labelsize=14)
        ax.set_xlabel('Time (day)', fontsize=14)
        ax.set_ylabel('BHP (bar)', fontsize=14)

        # Place the legend outside of the plot
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.grid()
        plt.tight_layout()  # Adjust the layout

        if savePlot:
            plt.savefig(f'{self.output_dir}/{fileName}.png', format='png', dpi=300)
        if showPlot:
            plt.show()

    def plot_pressureVsTime_semiLog(self, p_analytical=[], t_analytical=[], plotAnalytical=True,  showPlot=True, savePlot=False, fileName='P_logT'):
        # Plot the results
        fig, ax = plt.subplots(figsize=(11,8))

        # Loop over the results
        for (color, (keys, df)) in zip(self.colors, self.results.items()):
            # Plot data
            label = ', '.join(f'{label}: {value}' for label, value in zip(self.labels, keys))
            ax.plot(df['time'][:], df['PRD : BHP (bar)'][:], marker='o', linestyle='-', label=label, color=color)

        if plotAnalytical:
            ax.plot(t_analytical, p_analytical, color='r')

        ax.tick_params(labelsize=14)
        ax.set_xlabel('Time (day)', fontsize=14)
        ax.set_ylabel('BHP (bar)', fontsize=14)
        ax.set_xscale('log')

        # Place the legend outside of the plot
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.grid()
        plt.tight_layout()  # Adjust the layout
        if savePlot:
            plt.savefig(f'{self.output_dir}/{fileName}.png', format='png', dpi=300)
        if showPlot:
            plt.show()

    def plot_pD_tD(self, tD_analytical=[], pD_analytical=[], plotAnalytical=True, showPlot=True, savePlot=False,
                   fileName='PD_TD'):
        # Plot the results
        fig, ax = plt.subplots(figsize=(11, 8))

        # Loop over the results
        for (color, (keys, df)) in zip(self.colors, self.results.items()):
            # Plot data
            tD_simulation = DimensionlessTime(self.inputs, df['time'])
            pD_simulation = DimensionlessPressure(self.inputs, df['PRD : BHP (bar)'])
            # Create a label from the keys
            label = ', '.join(f'{label}: {value}' for label, value in zip(self.labels, keys))
            ax.plot(tD_simulation, pD_simulation, marker='o', linestyle='-', label=label, color=color)

        if plotAnalytical:
            ax.plot(tD_analytical, pD_analytical, color='r')

        ax.tick_params(labelsize=14)
        ax.set_xlabel('TD', fontsize=14)
        ax.set_ylabel('PD', fontsize=14)
        # Place the legend outside of the plot
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid()
        plt.tight_layout()  # Adjust the layout
        if savePlot:
            plt.savefig(f'{self.output_dir}/{fileName}.png', format='png', dpi=300)
        if showPlot:
            plt.show()

    def plot_pD_tD_log(self, tD_analytical=[], pD_analytical=[], plotAnalytical=True, showPlot=True, savePlot=False,
                   fileName='PD_TD_log'):
        # Plot the results
        fig, ax = plt.subplots(figsize=(11, 8))

        # Loop over the results
        for (color, (keys, df)) in zip(self.colors, self.results.items()):
            # Plot data
            tD_simulation = DimensionlessTime(self.inputs, df['time'])
            pD_simulation = DimensionlessPressure(self.inputs, df['PRD : BHP (bar)'])
            # Create a label from the keys
            label = ', '.join(f'{label}: {value}' for label, value in zip(self.labels, keys))
            ax.plot(tD_simulation, pD_simulation, marker='o', linestyle='-', label=label, color=color)

        if plotAnalytical:
            ax.plot(tD_analytical, pD_analytical, color='r')

        ax.tick_params(labelsize=14)
        ax.set_xlabel('TD', fontsize=14)
        ax.set_ylabel('PD', fontsize=14)
        ax.set_xscale('log')
        # Place the legend outside of the plot
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid()
        plt.tight_layout()  # Adjust the layout
        if savePlot:
            plt.savefig(f'{self.output_dir}/{fileName}.png', format='png', dpi=300)
        if showPlot:
            plt.show()

    def showVtk(self, fileName):
        # get vts data
        mesh = pv.read(f'{self.output_dir}/{fileName}.vts')

        # define plotter
        plotter = pv.Plotter()

        # set temperature as active scalar
        mesh.set_active_scalars('pressure')

        # add threshold levels
        thresT = mesh.threshold([0, self.inputs['Pi']], invert=False)

        # add outline of mesh
        outline = mesh.outline()

        # add elements to plotter
        plotter.set_background('#52576c')
        plotter.add_mesh(outline, color='k')
        plotter.add_mesh(thresT, cmap='coolwarm', show_edges=True,
                         scalar_bar_args={'title': 'pressure', 'fmt': '%.4f'})

        _ = plotter.add_axes(line_width=5, labels_off=False)

        plotter.show()

    def imageVtk(self, fileName, pressure_range=None,  showImage=False):
        # get vts data
        mesh = pv.read(f'{self.output_dir}\\{fileName}.vtk')

        # set temperature as active scalar
        mesh.set_active_scalars('pressure')

        # add threshold levels
        thresT = mesh.threshold([0, self.inputs['Pi']], invert=False)

        # add outline of mesh
        outline = mesh.outline()

        # define plotter
        plotter = pv.Plotter(off_screen=True)

        # add elements to plotter
        plotter.set_background('#52576c')
        plotter.add_mesh(outline, color='k')
        # If pressure range is specified, use it for color limits
        if pressure_range is not None:
            plotter.add_mesh(thresT, cmap='coolwarm', show_edges=True, clim=pressure_range,
                             scalar_bar_args={'title': 'pressure', 'fmt': '%.4f'})
        else:
            plotter.add_mesh(thresT, cmap='coolwarm', show_edges=True,
                             scalar_bar_args={'title': 'pressure', 'fmt': '%.4f'})

        _ = plotter.add_axes(line_width=5, labels_off=False)
        outputDir_vtk = f'{self.output_dir}/VTKs'
        if not os.path.exists(outputDir_vtk):
            os.makedirs(outputDir_vtk)

        # Render and create a screenshot
        plotter.screenshot(f'{outputDir_vtk}/{fileName}.png')

        if showImage:
            # Display the image
            Image(filename=f'{outputDir_vtk}/{fileName}.png')

    def makeGif(self, num_report_steps, output_file, duration):
        images = []
        outputDir_vtk = f'{self.output_dir}/VTKs'
        for ith_step in range(num_report_steps):
            filename = os.path.join(outputDir_vtk, 'solution' + str(ith_step) + '.png')
            img = PILImage.open(filename)
            images.append(img)
        images[0].save(f'{self.output_dir}/{output_file}', save_all=True, append_images=images[1:],
                        duration=duration , loop=0)
