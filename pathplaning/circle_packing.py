import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from shapely.geometry import Polygon, Point

# Function to calculate circle centers within a bounding box of the polygon
def calculate_circle_centers_polygon(vertices, r):
    centers = []

    # Create Shapely polygon object
    polygon = Polygon(vertices)

    # Get bounding box for the entire polygon
    x_coords, y_coords = zip(*vertices)
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)


    # Saftey factor one additional 1 is considered
    xw = abs(x_max - x_min) + 1
    yw = abs(y_max - y_min) + 1

    # Calculate grid dimensions for circle placement
    m = (int((xw - r) / (1.5 * r)) + 1) if ((xw - r) % (1.5 * r)) == (1.5 * r) else (int((xw - r) / (1.5 * r)) + 2)
    n1 = int((yw / (np.sqrt(3) * r)) - (np.sqrt(3) / 2)) + 2
    n2 = n1 if (yw / (np.sqrt(3) * r)) % 1 <= 0.5 else n1 - 1

    for l in range(1, m + 1):
        if l % 2 == 1:  # Odd column
            for k in range(1, n1 + 1):
                x = ((1.5 * l) - 1) * r + x_min
                y = round((k - 1) * np.sqrt(3) * r, 2) + y_min
                circle = Point(x, y ).buffer(r)  # Shapely circle as a polygon
                if circle.intersects(polygon):  # Check if circle intersects with polygon
                    centers.append((x, y))
        else:  # Even column
            for k in range(1, n2 + 1):
                x = ((1.5 * l) - 1) * r + x_min
                y = round((k - 1) * np.sqrt(3) * r + (np.sqrt(3) / 2) * r, 2) + y_min
                circle = Point(x, y).buffer(r)
                if circle.intersects(polygon):
                    centers.append((x, y))

    return centers

# Function to plot circles and polygon
def plot_circles_and_polygon(centers, vertices, r):
    fig, ax = plt.subplots()
    
    # Plot the polygon
    polygon = MplPolygon(vertices, closed=True, fill=None, edgecolor='b', label='Polygon')
    ax.add_patch(polygon)
    
    # Plot the circles
    for (cx, cy) in centers:
        circle = plt.Circle((cx, cy), r, color='r', fill=False)
        ax.add_patch(circle)
        ax.plot(cx, cy, 'ko')  # Mark the center of the circle

    ax.set_aspect('equal', 'box')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Circles with Intersection Area in the Polygon')
    plt.legend()
    plt.grid(True)
    plt.show()
