from shapely.geometry import Point, Polygon

def intersection_area(circle, quadrilateral):
    """
    Calculate the intersection area between a circle and a quadrilateral.
    """
    intersection = circle.intersection(quadrilateral)
    return intersection.area if not intersection.is_empty else 0

def is_circle_intersecting_quadrilateral(center, radius, quad_vertices):
    """
    Check if a circle intersects with a quadrilateral.
    """
    circle = Point(center).buffer(radius)  # Create a circle with given center and radius
    quadrilateral = Polygon(quad_vertices)  # Define the quadrilateral
    return circle.intersects(quadrilateral), circle

def get_valid_circle_centers(radius, centers, quad_vertices_list, target_index):
    """
    Return the centers of circles that intersect with the specified target quadrilateral,
    and append the center to the target quadrilateral if it has the maximum intersection area.
    """
    target_quad_vertices = quad_vertices_list[target_index]
    valid_centers = []
    
    for center in centers:
        # Check if circle intersects the target quadrilateral and get the circle object
        intersects, circle = is_circle_intersecting_quadrilateral(center, radius, target_quad_vertices)
        
        if intersects:
            max_intersection_area = 0
            max_quad_index = -1

            # Calculate intersection area with all quadrilaterals
            for i, quad_vertices in enumerate(quad_vertices_list):
                quadrilateral = Polygon(quad_vertices)
                intersection_area_value = intersection_area(circle, quadrilateral)
                
                # Track the maximum intersection area and its respective quadrilateral index
                if intersection_area_value > max_intersection_area:
                    max_intersection_area = intersection_area_value
                    max_quad_index = i
            
            # Append the center to the target quadrilateral if it has the maximum intersection area
            if max_quad_index == target_index:
                valid_centers.append(center)
                
    return valid_centers