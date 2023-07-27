# Import the gmsh library. Gmsh is a powerful mesh generation tool
# with a rich set of features for geometry creation and manipulation.
import gmsh

# Initialize the Gmsh environment. This is a necessary step that sets up
# the internal data structures Gmsh uses.
gmsh.initialize()

# Create a new model named "02Mesh". In Gmsh, a model is a container that
# can hold geometric entities (like points, lines, surfaces), as well as
# meshes, post-processing data, and more.
gmsh.model.add("02Mesh")

# Set the characteristic length (or mesh size). This is a parameter that
# controls the size of the elements in the mesh. A smaller value will result
# in a finer mesh, while a larger value will result in a coarser mesh.
lc = 0.1

# Define the corners of the square in the x-y plane. Each point is defined
# by its x, y, and z coordinates. The fourth argument to addPoint is the
# target mesh size near the point. This allows you to control the mesh size locally.
p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(1, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(1, 1, 0, lc)
p4 = gmsh.model.geo.addPoint(0, 1, 0, lc)

# Define the edges of the square. Each edge is a line connecting two points.
# These lines will form the boundaries of our surface.
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# Define a curve loop. A curve loop is a loop that is formed by connecting
# multiple curves end-to-end. In this case, the curve loop forms the boundary
# of the square. The order of the lines in the list is important as it determines
# the direction of the loop.
curve_loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])

# Define a plane surface bounded by the curve loop. This creates a 2D surface
# that can be meshed. In the context of geological modeling, this could represent
# a layer or stratum.
surface = gmsh.model.geo.addPlaneSurface([curve_loop])

# Extrude the surface in the z-direction to create a 2.5D mesh. The arguments
# to extrude specify the translation vector for the extrusion. This step gives
# thickness to our previously 2D surface, turning it into a 3D volume. However,
# since the meshing is still primarily in the 2D plane, we refer to this as a 2.5D mesh.
#dimTags: A list of tuples, where each tuple contains two integers. The first integer is
#the dimension of the entity (0 for points, 1 for lines, 2 for surfaces), and the second
# integer is the tag of the entity. For example, [(2, surface)] would extrude the surface
# with the given tag.
#dx, dy, dz: These are the translation values in the x, y, and z directions, respectively.
# These values determine the direction and length of the extrusion.

gmsh.model.geo.extrude([(2, surface)], 0, 0, 0.1)

# Synchronize the Gmsh model with the underlying geometric model. This is a crucial
# step that updates the Gmsh model with the changes we've made to the geometry.
# It must be done before generating the mesh.
gmsh.model.geo.synchronize()

# Generate the 3D mesh. This creates the actual mesh based on the geometry and
# settings we've defined. The argument '3' tellsGmsh to create a 3D mesh.
gmsh.model.mesh.generate(3)

# Write the mesh to a file named "02Mesh". This allows us to save our work and use
# the mesh in other programs or simulations.
gmsh.write("02Mesh.msh")

# Run the Gmsh GUI to visualize the mesh. This opens a new window where you can
# inspect the mesh. Visualization is a crucial part of the meshing process, as it
# allows you to check the quality of the mesh and identify any issues.
gmsh.fltk.run()

# Finalize the Gmsh environment. This should be done when you're finished using
# the gmsh library. It cleans up the internal data structures and frees the
# resources used by Gmsh.
gmsh.finalize()