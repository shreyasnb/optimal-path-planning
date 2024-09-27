import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from matplotlib.path import Path
import pandas as pd
# Load the data points from the CSV file
data_points= np.genfromtxt('points.csv', delimiter=',')
def reduce_to_polygon(hull_points, num_sides):
    current_sides = len(hull_points)
    if current_sides == num_sides:
        return hull_points  # Already has the correct number of points
    
    if current_sides > num_sides:
        # Reduce the number of sides by selecting evenly spaced points
        step = current_sides / num_sides
        reduced_points = [hull_points[int(i * step) % current_sides] for i in range(num_sides)]
        return np.array(reduced_points)
    
    # If fewer points, interpolate additional points
    return interpolate_points(hull_points, num_sides)

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

def scale_polygon(polygon_points, points, scale_factor=0.1):
    centroid = np.mean(polygon_points, axis=0)
    scaled_polygon = (polygon_points - centroid) * scale_factor + centroid
    
    # Increase scale factor until all points are inside the polygon
    while not all_inside_polygon(scaled_polygon, points):
        scale_factor += 0.05
        scaled_polygon = (polygon_points - centroid) * scale_factor + centroid
    
    return scaled_polygon

def all_inside_polygon(polygon, points):
    path = Path(polygon)
    return np.all(path.contains_points(points))

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


def ensure_10_sided_polygon(polygon_points, num_sides=10):
    current_sides = len(polygon_points)
    
    if current_sides == num_sides:
        return polygon_points
    
    if current_sides < num_sides:
        return interpolate_points(polygon_points, num_sides)
    
    if current_sides > num_sides:
        step = current_sides / num_sides
        reduced_points = [polygon_points[int(i * step) % current_sides] for i in range(num_sides)]
        return np.array(reduced_points)


def area_of_polygon(points):
    n = len(points)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += points[i][0] * points[j][1] - points[j][0] * points[i][1]
    New_area=abs(area)/2
    return New_area
# Function to create an array of points representing the sides of the polygon
def create_polygon_array(polygon_points):
    return np.concatenate((polygon_points, polygon_points[:1]))  # Close the polygon

# Function to find the closest data point for each vertex of the polygon
def find_closest_data_points(polygon_points, data_points):
    closest_points = []
    for point in polygon_points:
        distances = np.linalg.norm(data_points - point, axis=1)
        closest_points.append(data_points[np.argmin(distances)])
    return np.array(closest_points)

# Function to check if any convex hull vertex lies on any side of the polygon
def check_hull_on_polygon_sides(hull_points, polygon_points):
    active_indices = [True] * len(polygon_points)
    for hull_point in hull_points:
        for i in range(len(polygon_points) - 1):
            side_start = polygon_points[i]
            side_end = polygon_points[i + 1]
            if is_point_on_line_segment(hull_point, side_start, side_end):
                active_indices[i] = False
                active_indices[i + 1] = False
    return active_indices

# Helper function to determine if a point lies on a line segment
def is_point_on_line_segment(pt, line_start, line_end):
    cross_product = np.cross(line_end - line_start, pt - line_start)
    if abs(cross_product) > 1e-10:  # Close enough to zero
        return False
    dot_product = np.dot(pt - line_start, pt - line_end)
    return dot_product <= 0

# Function to check if any side of the polygon crosses a vertex of the convex hull
def check_if_crossed(hull_points, initial_polygon_points, final_polygon_points):
    crossed = False
    crossing_sides = []

    for i in range(len(initial_polygon_points) ):
        for hull_point in hull_points:
            # Get the vertices of the polygon side
            side_start_initial = initial_polygon_points[i]
            side_start_final = final_polygon_points[i]
            if i==len(initial_polygon_points)-1:
                side_end_initial = initial_polygon_points[ 0]
                side_end_final = final_polygon_points[0]
            else:
                side_end_initial = initial_polygon_points[i + 1]
                side_end_final = final_polygon_points[i + 1]

            # Calculate the vectors
            vector_side_initial = side_end_initial - side_start_initial
            vector_hull_initial = hull_point - side_start_initial
            vector_side_final = side_end_final - side_start_final
            vector_hull_final = hull_point - side_start_final
            
            # Compute the cross products
            cross_initial = np.cross(vector_side_initial, vector_hull_initial)
            cross_final=np.cross(vector_side_final, vector_hull_final)
            # Check if the signs of the cross products are different
            if cross_initial * cross_final <= 0:
                crossed = True
                crossing_sides.append(i)
                if i==len(initial_polygon_points)-1:
                    crossing_sides.append(0)
                else:
                    crossing_sides.append(i+1)
                break  # Stop checking this hull point for this side

    return crossed, crossing_sides


# Main movement logic function
def move_polygon_vertices(polygon_points, closest_points, hull_points, max_move_factor=1.0, step=0.009):
    move_factors = [0.0] * len(polygon_points)  # Initialize move factors
    active_indices = check_hull_on_polygon_sides(data_points, polygon_points)

    # Capture the initial polygon points for the first crossing check
    initial_polygon_points = polygon_points.copy()

    while True:
        # Step 1: Increase move factors for active vertices
        for i in range(len(polygon_points)):
            if active_indices[i]:  # Only increase for active indices
                move_factors[i] += step
                if move_factors[i] > max_move_factor :
                    move_factors[i] = max_move_factor
                    active_indices[i] = False
                    if i==len(polygon_points)-1:
                        active_indices[0]=False
        
                      # Cap at max move factor

        # Step 2: Move the polygon vertices towards the closest points
        for i in range(len(polygon_points)):
            if active_indices[i]:
                direction = closest_points[i] - polygon_points[i]
                if np.linalg.norm(direction) > 0:
                    # Move the point towards the closest point proportionally to the move factor
                    polygon_points[i] += (move_factors[i] / max_move_factor) * direction / np.linalg.norm(direction)

        # Step 3: Check if any side has crossed a convex hull vertex
        crossed, crossing_sides = check_if_crossed(data_points, initial_polygon_points, polygon_points)
        if crossed:
            for index in crossing_sides:
                active_indices[index] = False
                active_indices  # Mark these indices as inactive

        # Break the loop if all indices are inactive
        if all(not active for active in active_indices):
            break

    return polygon_points


# Main function to run the entire process
def main():
    hull = ConvexHull(data_points)
    hull_points = data_points[hull.vertices]

    # Reduce the hull to a 10-sided polygon
    polygon_points = reduce_to_polygon(hull_points, 10)

    # Enforce concavity if needed
    if len(polygon_points) < 10:
        concave_polygon_points = enforce_concavity(polygon_points)
        scaled_polygon_points = scale_polygon(concave_polygon_points, data_points)
        final_polygon_points = add_points_with_angles(scaled_polygon_points)


  
    else:
        # If the polygon already has 10 sides, just scale and plot
        scaled_polygon_points = scale_polygon(polygon_points, data_points)
        final_polygon_points = ensure_10_sided_polygon(scaled_polygon_points, 10)

    

    
    # Find closest points in the dataset for each vertex
    closest_points = find_closest_data_points(final_polygon_points, data_points)

    # Move vertices toward closest data points
    final_polygon_points = move_polygon_vertices(final_polygon_points, closest_points, hull_points)

    # Plotting results
    plt.figure(figsize=(8, 6))
    plt.plot(data_points[:, 0], data_points[:, 1], 'o', label='Data Points')
    plt.fill(hull_points[:, 0], hull_points[:, 1], 'b-', alpha=0.3, label='Convex Hull')
    plt.fill(final_polygon_points[:, 0], final_polygon_points[:, 1], 'r-', alpha=0.3, label='Optimized Polygon')
    plt.plot(final_polygon_points[:, 0], final_polygon_points[:, 1], 'r-', lw=2)
    plt.scatter(final_polygon_points[:, 0], final_polygon_points[:, 1], c='red', label='Polygon Vertices')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Optimized Polygon and Convex Hull')
    plt.legend()
    plt.grid(True)
    plt.show()
    print("Final Polygon Points:", final_polygon_points)

if __name__ == "__main__":
    main()

