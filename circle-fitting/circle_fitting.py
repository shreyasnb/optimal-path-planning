import numpy as np
import matplotlib.pyplot as plt

# Function to calculate circle centers based on the formula provided
def calculate_circle_centers(x1, y1, x2, y2, r):
    # Calculate the width and height of the rectangle
    xw = abs(x2 - x1)
    yw = abs(y2 - y1)

    # Determine the bottom-left corner for reference
    x_min = min(x1, x2)
    y_min = min(y1, y2)

    # Calculate m and n1, n2
    m = (int((xw - r) / (1.5 * r)) + 1) if (xw - r) % (1.5 * r) == 1.5 * r else (int((xw - r) / (1.5 * r)) + 2)
    n1 = int((yw / (np.sqrt(3) * r)) - (np.sqrt(3) / 2)) + 2
    n2 = n1 if (yw / (np.sqrt(3) * r)) % 1 <= 0.5 else n1 - 1
    print(m,n1,n2)

    centers = []
    
    # Iterate through rows and columns to calculate circle centers
    for l in range(1, m + 1):

        if xw/yw == 2 or yw/xw == 2:
            if l % 2 == 1:  # Odd column

                for k in range(1, n1 + 2):
            
                    x = ((1.5 * l) - 1) * r
                    y = round((k - 1) * np.sqrt(3) * r,2)
                    centers.append((x, y))
                    
            else:  # Even column

                for k in range(1, n2 + 2):    
                    x = ((1.5 * l) - 1) * r
                    y = round((k - 1) * np.sqrt(3) * r + (np.sqrt(3) / 2) * r,2)
                    centers.append((x, y))
            # if x <= xw and y <= yw:  # Ensure the circle is within the rectangle
            # centers.append((x, y))

        else:
            if l % 2 == 1:  # Odd column

                for k in range(1, n1 + 1):
            
                    x = ((1.5 * l) - 1) * r
                    y = round((k - 1) * np.sqrt(3) * r,2)
                    centers.append((x, y))
            else:  # Even column

                for k in range(1, n2 + 1):    
                    x = ((1.5 * l) - 1) * r
                    y = round((k - 1) * np.sqrt(3) * r + (np.sqrt(3) / 2) * r,2)
                    centers.append((x, y))
            # if x <= xw and y <= yw:  # Ensure the circle is within the rectangle
            # centers.append((x, y))
    print(centers)
    return centers

# Function to plot the circles and rectangle
def plot_circles_and_rectangle(centers, x1, y1, x2, y2, r):
    fig, ax = plt.subplots()
    plt.figure
    # Plot the rectangle
    rect_x = [x1, x2, x2, x1, x1]
    rect_y = [y1, y1, y2, y2, y1]
    ax.plot(rect_x, rect_y, 'b-', label='Rectangle')
    
    # Plot the circles
    for (cx, cy) in centers:
        circle = plt.Circle((cx, cy), r, color='r', fill=False)
        ax.add_patch(circle)
        ax.plot(cx, cy, 'ko')  # Mark the center of the circle

    ax.set_aspect('equal', 'box')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Circles inside the Rectangle')
    plt.legend()
    plt.grid(True)
    plt.show()

# Main function to execute the process
def main():
    # Input rectangle diagonal coordinates
    x1, y1 = map(float, input("Enter the first diagonal coordinate of the rectangle (x1 y1): ").split())
    x2, y2 = map(float, input("Enter the second diagonal coordinate of the rectangle (x2 y2): ").split())
    r = float(input("Enter the radius of the circles: "))

    # Get circle centers
    centers = calculate_circle_centers(x1, y1, x2, y2, r)
    
    # Plot the rectangle and circles
    plot_circles_and_rectangle(centers, x1, y1, x2, y2, r)

# Execute the main function
if __name__ == "__main__":
    main()

