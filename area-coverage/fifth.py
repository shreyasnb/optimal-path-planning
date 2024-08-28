import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.optimize import minimize

# Define the 6 initial points
initial_points = np.array([
    [1, 1], [3, 1], [4, 3], [3, 5], [1, 5], [0, 3]
])

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

def create_concave_polygon(points, num_sides):
    # Generate interpolated points
    polygon_points = interpolate_points(points, num_sides)
    
    # Random perturbations to ensure points are not collinear
    perturbation_strength = 0.1
    perturbed_points = polygon_points + np.random.uniform(-perturbation_strength, perturbation_strength, polygon_points.shape)
    
    # Optimization to induce concavity
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
    
    return final_polygon_points

# Generate a 10-sided concave polygon
concave_polygon_points = create_concave_polygon(initial_points, 10)

# Plot the initial points and the resulting concave 10-sided polygon
plt.figure(figsize=(8, 6))
plt.plot(initial_points[:, 0], initial_points[:, 1], 'o', label='Initial Points')
plt.plot(concave_polygon_points[:, 0], concave_polygon_points[:, 1], 'r-', lw=2, label='Concave 10-sided Polygon')
plt.fill(concave_polygon_points[:, 0], concave_polygon_points[:, 1], 'r-', alpha=0.3)
plt.scatter(concave_polygon_points[:, 0], concave_polygon_points[:, 1], c='red', label='Polygon Vertices')

plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Concave 10-sided Polygon Generated from 6 Points')
plt.legend()
plt.grid(True)
plt.show()
