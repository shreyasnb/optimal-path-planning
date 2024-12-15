import matplotlib.pyplot as plt
from collections import defaultdict


def zigzag_raster_scan(nodes):
    # Step 1: Sort nodes by x first, then y to prioritize vertical sorting
    nodes.sort(key=lambda node: (node[0], node[1]))

    # Step 2: Group nodes by their x-coordinate (vertical rows)
    columns = defaultdict(list)
    for x, y in nodes:
        columns[x].append((x, y))

    # Step 3: Create the zig-zag path for vertical scanning
    zigzag_path = []
    for i, (x, col_nodes) in enumerate(sorted(columns.items())):
        # Reverse every other column to create zig-zag vertically
        if i % 2 == 1:
            col_nodes = col_nodes[::-1]
        zigzag_path.extend(col_nodes)

    return zigzag_path
