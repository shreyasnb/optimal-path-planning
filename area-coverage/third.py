import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

points=np.genfromtxt('pacman.csv', delimiter=',')

# Calculate the convex hull
hull = ConvexHull(points)
hull_points = points[hull.vertices]

# Function to reduce/simplify the convex hull to a 10-sided polygon
def reduce_to_polygon(hull_points, num_sides):
    current_sides = len(hull_points)
    if current_sides == num_sides:
        return hull_points  # Already has the correct number of points
    
    if current_sides > num_sides:
        # Reduce the number of sides by selecting evenly spaced points
        step = current_sides / num_sides
        reduced_points = [hull_points[int(i * step) % current_sides] for i in range(num_sides)]
        return np.array(reduced_points)
    
    # If fewer points, interpolate additional points (previous method)
    return interpolate_points(hull_points, num_sides)

# Function to interpolate additional points if fewer than required
def interpolate_points(hull_points, num_sides):
    current_sides = len(hull_points)
    new_points = []
    for i in range(current_sides):
        new_points.append(hull_points[i])
        next_point = hull_points[(i + 1) % current_sides]
        num_interpolations = (num_sides - current_sides) // current_sides
        for j in range(1, num_interpolations + 1):
            interpolated_point = hull_points[i] + (next_point - hull_points[i]) * (j / (num_interpolations + 1))
            new_points.append(interpolated_point)
        if len(new_points) >= num_sides:
            break
    return np.array(new_points[:num_sides])

# Function to scale the polygon outwards
def scale_polygon(polygon_points, points, scale_factor=1.1):
    centroid = np.mean(polygon_points, axis=0)
    scaled_polygon = (polygon_points - centroid) * scale_factor + centroid
    
    # Increase scale factor until all points are inside the polygon
    while not all_inside_polygon(scaled_polygon, points):
        scale_factor += 0.05
        scaled_polygon = (polygon_points - centroid) * scale_factor + centroid
    
    return scaled_polygon

# Check if all points are inside the polygon using cross products
def all_inside_polygon(polygon, points):
    from matplotlib.path import Path
    path = Path(polygon)
    return np.all(path.contains_points(points))

# Reduce or expand to form a 10-sided polygon
polygon_points = reduce_to_polygon(hull_points, 10)

# Scale the polygon to ensure all points are inside
scaled_polygon_points = scale_polygon(polygon_points, points)

# Plotting the data points and the scaled 10-sided polygon
plt.figure(figsize=(8, 6))
plt.plot(points[:, 0], points[:, 1], 'o', label='Data Points')
plt.fill(scaled_polygon_points[:, 0], scaled_polygon_points[:, 1], 'r-', alpha=0.3, label='Scaled 10-sided Polygon')
plt.plot(scaled_polygon_points[:, 0], scaled_polygon_points[:, 1], 'r-', lw=2)
plt.scatter(scaled_polygon_points[:, 0], scaled_polygon_points[:, 1], c='red', label='Polygon Vertices')

plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Scaled 10-sided Polygon Enclosing All Points')
plt.legend()
plt.show()
