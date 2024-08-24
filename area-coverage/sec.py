import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

# Example n-sided polygon vertices
polygon_vertices = np.array([
    (1, 0), (0.995, 0.099), (0.98, 0.197), (0.955, 0.294),
    (0.921, 0.388), (0.877, 0.477), (0.824, 0.56), (0.762, 0.637),
    (0.691, 0.706), (0.612, 0.767)
])

# Compute the convex hull
hull = ConvexHull(polygon_vertices)
hull_points = polygon_vertices[hull.vertices]

# Compute the diameter of the smallest circle enclosing the convex hull
def max_distance(points):
    max_dist = 0
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            dist = np.linalg.norm(points[i] - points[j])
            if dist > max_dist:
                max_dist = dist
    return max_dist

diameter = max_distance(hull_points)
radius = diameter / 2

# Generate decagon vertices
num_sides = 10
angles = np.linspace(0, 2 * np.pi, num_sides, endpoint=False)
decagon_vertices = np.array([(radius * np.cos(angle), radius * np.sin(angle)) for angle in angles])

# Plot the results
plt.figure(figsize=(10, 10))
plt.plot(*zip(*hull_points, hull_points[0]), 'g-', label='Convex Hull')
plt.fill(*zip(*hull_points, hull_points[0]), 'g', alpha=0.3)
plt.plot(*zip(*decagon_vertices, decagon_vertices[0]), 'b-', label='10-Sided Polygon')
plt.scatter(*zip(*polygon_vertices), c='r', s=10, label='Original Points')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Decagon Encompassing n-Sided Polygon')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()