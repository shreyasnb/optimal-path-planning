import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.optimize import minimize

def interpolate_points(points, num_sides):
    current_sides = len(points)
    if current_sides >= num_sides:
        raise ValueError("Number of initial points must be less than the desired number of sides.")
    
    interpolated_points = []
    new_points_needed = num_sides - current_sides
    
    # Generate interpolated points
    for i in range(current_sides):
        interpolated_points.append(points[i])
        next_point_idx = (i + 1) % current_sides
        segment_length = (new_points_needed + current_sides - 1) // current_sides
        for j in range(1, segment_length + 1):
            ratio = j / (segment_length + 1)
            new_point = points[i] + ratio * (points[next_point_idx] - points[i])
            interpolated_points.append(new_point)
            new_points_needed -= 1
            if new_points_needed <= 0:
                break
    
    return np.array(interpolated_points)

def angle_between(v1, v2):
    """Calculate the angle between two vectors in degrees."""
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    return np.degrees(angle)

def enforce_concavity(polygon_points, concave_angle):
    """Adjusts points to enforce concavity at certain vertices."""
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
            polygon_points[i] -= 0.5 * adjustment_vector  # Move the point inward to increase the concavity
            
    return polygon_points

def scale_polygon(polygon_points, scaling_factor):
    """Scale the polygon outward from its centroid."""
    centroid = np.mean(polygon_points, axis=0)
    scaled_points = centroid + scaling_factor * (polygon_points - centroid)
    return scaled_points

def is_point_inside_polygon(point, polygon):
    """Check if a point is inside a polygon using the ray-casting algorithm."""
    x, y = point
    n = len(polygon)
    inside = False
    p1x, p1y = polygon[0]
    for i in range(n + 1):
        p2x, p2y = polygon[i % n]
        if min(p1y, p2y) < y <= max(p1y, p2y):
            if x <= max(p1x, p2x):
                if p1y != p2y:
                    xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                if p1x == p2x or x <= xinters:
                    inside = not inside
        p1x, p1y = p2x, p2y
    return inside

def expand_polygon_to_include_points(polygon_points, initial_points):
    """Expand the polygon until all initial points are inside the polygon."""
    scaling_factor = 1.0
    while not all(is_point_inside_polygon(point, polygon_points) for point in initial_points):
        scaling_factor += 0.1
        polygon_points = scale_polygon(polygon_points, scaling_factor)
    return polygon_points

def create_concave_polygon(initial_points, num_sides, concave_angle=150):
    # Generate interpolated points
    polygon_points = interpolate_points(initial_points, num_sides)
    
    # Random perturbations to ensure points are not collinear
    perturbation_strength = 0.1
    perturbed_points = polygon_points + np.random.uniform(-perturbation_strength, perturbation_strength, polygon_points.shape)
    
    # Enforce concavity
    perturbed_points = enforce_concavity(perturbed_points, concave_angle)
    
    # Optimization to fine-tune concavity and area
    def objective(x):
        poly_points = x.reshape((-1, 2))
        hull = ConvexHull(poly_points)
        area = hull.volume
        return area

    def constraint(x):
        poly_points = x.reshape((-1, 2))
        hull = ConvexHull(poly_points)
        hull_points = poly_points[hull.vertices]
        return len(hull_points) - num_sides
    
    initial_guess = perturbed_points.flatten()
    constraints = {'type': 'eq', 'fun': constraint}
    bounds = [(None, None)] * len(initial_guess)
    
    result = minimize(objective, initial_guess, method='SLSQP', constraints=constraints, bounds=bounds)
    final_polygon_points = result.x.reshape((-1, 2))
    
    # Expand the polygon to ensure all initial points lie inside
    final_polygon_points = expand_polygon_to_include_points(final_polygon_points, initial_points)
    
    return final_polygon_points

# Define the initial points (example)
initial_points = np.array([
    [1, 1], [3, 1], [4, 3], [3, 5], [1, 5], [0, 3]
])

# Specify the number of sides for the polygon
num_sides = 10  # Example: Generate a 10-sided polygon

# Generate the concave polygon
concave_polygon_points = create_concave_polygon(initial_points, num_sides)

# Plot the initial points and the resulting concave n-sided polygon
plt.figure(figsize=(8, 6))
plt.plot(initial_points[:, 0], initial_points[:, 1], 'o', label='Initial Points')
plt.plot(concave_polygon_points[:, 0], concave_polygon_points[:, 1], 'r-', lw=2, label=f'Concave {num_sides}-sided Polygon')
plt.fill(concave_polygon_points[:, 0], concave_polygon_points[:, 1], 'r-', alpha=0.3)
plt.scatter(concave_polygon_points[:, 0], concave_polygon_points[:, 1], c='red', label='Polygon Vertices')

plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title(f'Concave {num_sides}-sided Polygon Generated from Initial Points')
plt.legend()
plt.grid(True)
plt.show()
