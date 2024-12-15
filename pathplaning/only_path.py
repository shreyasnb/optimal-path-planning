import numpy as np
import matplotlib.pyplot as plt
import circle_packing as cp
import center_filter as cf
import polygon as poly
import pathplan as pp

def get_paths_only(r):
    """Generate and return zig-zag paths without plotting."""
    
    # Generate quadrilaterals and circle centers
    quadrilaterals, common_vertex, vertices = poly.gen_quad()
    centers = cp.calculate_circle_centers_polygon(vertices, r)

    # Populate nodes inside each quadrilateral
    n = []
    for target_index, quad in enumerate(quadrilaterals):
        n.append(cf.get_valid_circle_centers(r, centers, quadrilaterals, target_index))

    # Generate zig-zag paths for each quadrilateral's nodes
    zigzag_paths = [pp.zigzag_raster_scan(nodes) for nodes in n]
    
    return zigzag_paths  # Return paths for external use

# Main code to get the paths without plotting
r = 1
R = 20
paths = get_paths_only(r)

# Display or use paths as needed
for i, path in enumerate(paths):
    print(f"Path {i + 1}: {path}")