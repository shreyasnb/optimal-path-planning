import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
points=np.genfromtxt('pacman.csv',delimiter=',')

# Calculate the convex hull
hull = ConvexHull(points)
hull_points = points[hull.vertices]

# Function to reduce the convex hull to a 10-sided polygon
def reduce_to_polygon(hull_points, num_sides):
    current_sides = len(hull_points)
    
    # If there are fewer points, interpolate between them to reach the desired number of sides
    if current_sides < num_sides:
        return interpolate_points(hull_points, num_sides)
    
    # Reduce the number of sides by selecting evenly spaced points
    step = current_sides / num_sides
    reduced_points = [hull_points[int(i * step) % current_sides] for i in range(num_sides)]
    return np.array(reduced_points)

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

# Function to enforce concavity by adjusting vertices
def enforce_concavity(polygon_points, concave_angle=150):
    def angle_between(v1, v2):
        """Calculate the angle between two vectors in degrees."""
        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        angle = np.arccos(np.clip(cos_theta, -1.0, 1.0))
        return np.degrees(angle)

    for i in range(len(polygon_points)):
        prev_point = polygon_points[i - 1]
        current_point = polygon_points[i]
        next_point = polygon_points[(i + 1) % len(polygon_points)]
        
        v1 = current_point - prev_point
        v2 = next_point - current_point
        angle = angle_between(v1, v2)
        
        # If the angle is less than the desired concave angle, adjust the points
        if angle < concave_angle:
            adjustment_vector = 0.5 * (v1 + v2)
            polygon_points[i] -= 0.5 * adjustment_vector  # Move the point inward to increase concavity
    
    return polygon_points

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

def add_points_with_angles(vertices, target_sides=10, offset_factor=-0.01):
    n = len(vertices)
    if n >= target_sides:
        raise ValueError("The polygon already has 10 or more sides.")
    
    points_to_add = target_sides - n
    new_vertices = []

    add_per_side = points_to_add // n
    extra_points = points_to_add % n
    
    for i in range(n):
        new_vertices.append(vertices[i])
        next_index = (i + 1) % n
        for j in range(add_per_side + (1 if extra_points > 0 else 0)):
            t = (j + 1) / (add_per_side + 1 + (1 if extra_points > 0 else 0))
            new_point = (1 - t) * np.array(vertices[i]) + t * np.array(vertices[next_index])
            
            # Offset the new point perpendicularly to create angles
            side_vector = np.array(vertices[next_index]) - np.array(vertices[i])
            perpendicular_vector = np.array([-side_vector[1], side_vector[0]])
            perpendicular_vector = perpendicular_vector / np.linalg.norm(perpendicular_vector)
            offset = offset_factor * np.linalg.norm(side_vector) * perpendicular_vector
            
            new_point_with_offset = new_point + offset
            new_vertices.append(new_point_with_offset.tolist())
        
        if extra_points > 0:
            extra_points -= 1

    return np.array(new_vertices)

# Reduce or expand to form a 10-sided polygon
polygon_points = reduce_to_polygon(hull_points, 10)

# Enforce concavity on the polygon
concave_polygon_points = enforce_concavity(polygon_points)

# Scale the polygon to ensure all points are inside
scaled_polygon_points = add_points_with_angles(scale_polygon(concave_polygon_points, points))


# Print the coordinates of the 10-sided polygon
print("Coordinates of the 10-sided polygon:")
print(scaled_polygon_points)

# Plotting the data points and the scaled 10-sided concave polygon
plt.figure(figsize=(8, 6))
plt.plot(points[:, 0], points[:, 1], 'o', label='Data Points')
plt.fill(scaled_polygon_points[:, 0], scaled_polygon_points[:, 1], 'r-', alpha=0.3, label='Scaled 10-sided Concave Polygon')
plt.plot(scaled_polygon_points[:, 0], scaled_polygon_points[:, 1], 'r-', lw=2)
plt.scatter(scaled_polygon_points[:, 0], scaled_polygon_points[:, 1], c='red', label='Polygon Vertices')


plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Scaled 10-sided Concave Polygon Enclosing All Points')
plt.legend()
plt.grid(True)
plt.show()