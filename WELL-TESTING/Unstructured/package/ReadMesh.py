import meshio
import matplotlib.pyplot as plt

# Read the .msh file
import gmsh

# Initialize the Gmsh environment
gmsh.initialize()

# Load the .msh file
gmsh.open('C:\meshes\mesh_60_real_6.msh')

# Run the Gmsh GUI to visualize the mesh
gmsh.fltk.run()

# Finalize the Gmsh environment
gmsh.finalize()
