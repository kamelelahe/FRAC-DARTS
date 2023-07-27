import matplotlib.pyplot as plt
import numpy as np
from ..equations.DimensionLess import DimensionlessPressure, DimensionlessTime
import matplotlib.cm as cm
import pyvista as pv
from IPython.display import Image

class Visualize:

    def __init__(self, results, inputs, labels, output_dir, MuList):
        self.output_dir = output_dir
        self.results = results
        self.inputs = inputs
        self.labels = labels
        self.MuList= MuList
        # Create a colormap
        self.colors = cm.rainbow(np.linspace(0, 1, len(self.results)))
        self.colorsAnalytical = cm.rainbow(np.linspace(0, 1, len(self.MuList)))

    def plot_pressureVsTime(self, p_analytical=[], t_analytical=[] ,plotAnalytical=True,  showPlot=True, savePlot=False, fileName='P_T'):
        # Plot the results
        fig, ax = plt.subplots(figsize=(11,8))
        # Loop over the results
        for (color, (keys, df)) in zip(self.colors, self.results.items()):
            # Plot data
            label = ', '.join(f'{label}: {value}' for label, value in zip(self.labels, keys))
            ax.plot(df['time'], df['PRD : BHP (bar)'], marker='o', linestyle='-', label=label, color=color)

        if plotAnalytical:
            for (color, (key, value)) in zip(self.colorsAnalytical, p_analytical.items()):
                label = 'Mu : ' + str(key)
                ax.plot(t_analytical, value, label=label, color=color)

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

    def plot_pressureVsTime_semiLog(self, p_analytical=[], t_analytical=[],  plotAnalytical=True,  showPlot=True, savePlot=False, fileName='P_logT'):
        # Plot the results
        fig, ax = plt.subplots(figsize=(11,8))

        # Loop over the results
        for (color, (keys, df)) in zip(self.colors, self.results.items()):
            # Plot data
            label = ', '.join(f'{label}: {value}' for label, value in zip(self.labels, keys))
            ax.plot(df['time'][:], df['PRD : BHP (bar)'][:], marker='o', linestyle='-', label=label, color=color)

        if plotAnalytical:
            for (color, (key, value)) in zip(self.colorsAnalytical, p_analytical.items()):
                label = 'Mu : ' + str(key)
                ax.plot(t_analytical, value, label=label, color=color)

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

    def plot_pD_tD(self, tD_analytical={}, pD_analytical={}, plotAnalytical=True, showPlot=True, savePlot=False,
                   fileName='PD_TD'):
        # Plot the results
        fig, ax = plt.subplots(figsize=(11, 8))
        # Loop over the results
        for (color, (keys, df)) in zip(self.colors, self.results.items()):
            # Plot data
            if 'Mu' in self.labels:
                mu_index = self.labels.index('Mu')
                self.inputs['mu'] = keys[mu_index]
            tD_simulation = DimensionlessTime(self.inputs, df['time'])
            pD_simulation = DimensionlessPressure(self.inputs, df['PRD : BHP (bar)'])
            # Create a label from the keys
            label = ', '.join(f'{label}: {value}' for label, value in zip(self.labels, keys))
            ax.plot(tD_simulation, pD_simulation, marker='o', linestyle='-', label=label, color=color)

        if plotAnalytical:
            for (color, key) in zip(self.colorsAnalytical, pD_analytical.keys()):
                label = 'Mu : ' + str(key)
                ax.plot(tD_analytical[key], pD_analytical[key], label=label, color=color)

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

    def plot_pD_tD_log(self, tD_analytical=[], pD_analytical={}, plotAnalytical=True, showPlot=True, savePlot=False,
                   fileName='PD_TD_log'):
        # Plot the results
        fig, ax = plt.subplots(figsize=(11, 8))

        # Loop over the results
        for (color, (keys, df)) in zip(self.colors, self.results.items()):
            # Plot data
            if 'Mu' in self.labels:
                mu_index = self.labels.index('Mu')
                self.inputs['mu'] = keys[mu_index]
            tD_simulation = DimensionlessTime(self.inputs, df['time'])
            pD_simulation = DimensionlessPressure(self.inputs, df['PRD : BHP (bar)'])
            # Create a label from the keys
            label = ', '.join(f'{label}: {value}' for label, value in zip(self.labels, keys))
            ax.plot(tD_simulation, pD_simulation, marker='o', linestyle='-', label=label, color=color)

        if plotAnalytical:
            for (color, key) in zip(self.colorsAnalytical, pD_analytical.keys()):
                label = 'Mu : ' + str(key)
                ax.plot(tD_analytical[key], pD_analytical[key], label=label, color=color)


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
        mesh = pv.read(f'vtk_data/{fileName}.vts')

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

    def imageVtk(self, fileName, output_dir, showImage=True):
        # Specify the file
        i = 0  # adjust this according to the specific solution file you want

        # get vts data
        mesh = pv.read(f'vtk_data\\{fileName}.vts')

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
        plotter.add_mesh(thresT, cmap='coolwarm', show_edges=True,
                         scalar_bar_args={'title': 'pressure', 'fmt': '%.4f'})


        _ = plotter.add_axes(line_width=5, labels_off=False)

        # Render and create a screenshot
        plotter.screenshot(f'{output_dir}/{fileName}.png')

        if showImage:
            # Display the image
            Image(filename=f'{output_dir}/{fileName}.png')