#!/usr/bin/env python3
"""
This program converts a GDSII 2D layout file to multiple 3D STL files that can
be visualized in an external program (e.g., Blender).

USAGE:
    - edit the "layerstack" variable in the "CONFIGURATION" section below
    - run "gdsiistl file.gds"
OUTPUT:
    - the files file.gds_layername1.stl, file.gds_layername2.stl, ...

The program takes one argument, a path to a GDSII file. It reads shapes from
each layer of the GDSII file, converts them to polygon boundaries, then makes
a triangle mesh for each GDSII layer by extruding the polygons to given sizes.
"""

import sys # read command-line arguments
import gdspy # open gds file
from stl import mesh # write stl file (python package name is "numpy-stl")
import numpy as np # fast math on lots of points
import triangle # triangulate polygons

# get the input file name
if len(sys.argv) < 2: # sys.argv[0] is the name of the program
    print("Error: need exactly one file as a command line argument.")
    sys.exit(0)
gdsii_file_path = sys.argv[1]

########## CONFIGURATION (EDIT THIS PART) #####################################

# choose which GDSII layers to use
layerstack = {
    # layernumber: (zmin, zmax, 'layername'),
    1: (0, 550, 'substrate'),
    3: (552, 592, 'soi'),
    6: (592, 593, 'metal'),
}

########## INPUT ##############################################################

print('reading GDSII file {}...'.format(gdsii_file_path))
gdsii = gdspy.GdsLibrary()
gdsii.read_gds(gdsii_file_path, units='import')
# see https://gdspy.readthedocs.io/en/stable/index.html for documentation

print('extracting polygons...')
layers = {} # array to hold all geometry, sorted into layers

cells = gdsii.top_level() # get all cells that aren't referenced by another
for cell in cells: # loop through cells to read paths and polygons

    # $$$CONTEXT_INFO$$$ is a separate, non-standard compliant cell added
    # optionally by KLayout to store extra information not needed here.
    # see https://www.klayout.de/forum/discussion/1026/very-
    # important-gds-exported-from-k-layout-not-working-on-cadence-at-foundry
    if cell.name == '$$$CONTEXT_INFO$$$':
        continue # skip this cell

    # combine will all referenced cells (instances, SREFs, AREFs, etc.)
    cell = cell.flatten()

    # loop through paths in cell
    for path in cell.paths:
        lnum = path.layers[0] # GDSII layer number
        # create empty array to hold layer polygons if it doesn't yet exist
        layers[lnum] = [] if not lnum in layers else layers[lnum]
        # add paths (converted to polygons) that layer
        for poly in path.get_polygons():
            layers[lnum].append((poly, None, False))

    # loop through polygons (and boxes) in cell
    for polygon in cell.polygons:
        lnum = polygon.layers[0] # same as before...
        layers[lnum] = [] if not lnum in layers else layers[lnum]
        for poly in polygon.polygons:
            layers[lnum].append((poly, None, False))

"""
At this point, "layers" is a Python dictionary structured as follows:

layers = {
   0 : [ ([[x1, y1], [x2, y2], ...], None, False), ... ]
   1 : [ ... ]
   2 : [ ... ]
   ...
}

Each dictionary key is a GDSII layer number (0-255), and the value of the
dictionary at that key (if it exists; keys were only created for layers with
geometry) is a list of polygons in that GDSII layer. Each polygon is a 2-tuple
whose first element is a list of points (2-element lists with x and y
coordinates), second element is None (for the moment; this will be used later),
and third element is False (whether the polygon is clockwise; will be updated).
"""

########## TRIANGULATION ######################################################

print('triangulating polygons...')

num_triangles = {} # will store the number of triangles for each layer

# loop through all layers
for layer_number, polygons in layers.items():

    # but skip layer if it won't be exported
    if not layer_number in layerstack.keys():
        continue

    num_triangles[layer_number] = 0

    # loop through polygons in layer
    for index, (polygon, _, _) in enumerate(polygons):

        num_polygon_points = len(polygon)

        # determine whether polygon points are CW or CCW
        area = 0
        for i, v1 in enumerate(polygon): # loop through vertices
            v2 = polygon[(i+1) % num_polygon_points]
            area += (v2[0]-v1[0])*(v2[1]+v1[1]) # integrate area
        clockwise = area > 0

        # GDSII implements holes in polygons by making the polygon edge
        # wrap into the hole and back out along the same line. However,
        # this confuses the triangulation library, which fills the holes
        # with extra triangles. Avoid this by moving each edge back a
        # very small amount so that no two edges of the same polygon overlap.
        delta = 0.001 # inset each vertex by this much
        for i, v1 in enumerate(polygon):
            # calculate unit normal 2D vectors to each edge of a vertex
            v2 = polygon[(i+1) % num_polygon_points]
            v3 = polygon[(i-1) % num_polygon_points]
            normal2 = np.array([v2[1]-v1[1], v1[0]-v2[0]])
            normal2 /= np.linalg.norm(normal2)
            normal3 = np.array([v1[1]-v3[1], v3[0]-v1[0]])
            normal3 /= np.linalg.norm(normal3)
            if clockwise:
                normal2 = -normal2
                normal3 = -normal3
            # move vertex back along both normals
            polygon[i] = [v1[0] - normal2[0]*delta - normal3[0]*delta,
                          v1[1] - normal2[1]*delta - normal3[1]*delta]

        # triangulate: compute triangles to fill polygon
        point_array = np.arange(num_polygon_points)
        edges = np.transpose(np.stack((point_array, np.roll(point_array, 1))))
        triangles = triangle.triangulate(dict(vertices=polygon,
                                              segments=edges), opts='p')
        if not 'triangles' in triangles.keys():
            triangles['triangles'] = []

        # each line segment will make two triangles (for a rectangle), and the polygon
        # triangulation will be copied on the top and bottom of the layer.
        num_triangles[layer_number] += num_polygon_points*2 + \
                                       len(triangles['triangles'])*2
        polygons[index] = (polygon, triangles, clockwise)

########## EXTRUSION ##########################################################

print('extruding polygons...')

# loop through all layers
for layer in layers:

    # but skip layer if it won't be exported
    if not layer in layerstack.keys():
        continue

    # make a list of triangles
    layer_mesh_data = np.zeros(num_triangles[layer], dtype=mesh.Mesh.dtype)
    layer_pointer = 0

    for index, (polygon, triangles, clockwise) in enumerate(layers[layer]):

        zmin, zmax, layername = layerstack[layer]

        points_i = polygon
        points_i_min = np.insert(points_i, 2, zmin, axis=1)
        points_i_max = np.insert(points_i, 2, zmax, axis=1)
        points_j_min = np.roll(points_i_min, 1, axis=0)
        points_j_max = np.roll(points_i_max, 1, axis=0)
        lefts = np.stack((points_i_min, points_j_min, points_i_max), axis=1)
        rights = np.stack((points_j_min, points_j_max, points_i_max), axis=1)
        vs = triangles['vertices']
        ts = triangles['triangles']

        face_tris = np.take(vs, ts, axis=0)
        top = np.insert(face_tris, 2, zmax, axis=2)
        bottom = np.insert(face_tris, 2, zmin, axis=2)

        faces = np.concatenate((lefts, rights, top, bottom), axis=0)
        #faces = np.concatenate((lefts, rights), axis=0)
        num = len(faces)
        layer_mesh_data['vectors'][layer_pointer:(layer_pointer+len(faces))] = faces
        layer_pointer += len(faces)

    # save layer to STL file
    print('writing layer {} to file...'.format(layername))
    layer_mesh_object = mesh.Mesh(layer_mesh_data, remove_empty_areas=False)
    layer_mesh_object.save(gdsii_file_path + '_{}.stl'.format(layername))
