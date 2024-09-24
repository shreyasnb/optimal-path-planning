# Circle Fitting Algorithm in a Rectangle

This project implements a Python script to calculate the centers of circles packed inside a rectangle, given the rectangle's diagonal coordinates and the radius of the circles. The circles are plotted using `matplotlib`, and the algorithm ensures that the circles fit inside the rectangle.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Function Descriptions](#function-descriptions)
- [Example](#example)

## Overview
The script calculates the optimal positions of circles packed inside a rectangle using hexagonal packing, one of the most space-efficient methods. The circle radius and rectangle dimensions are given as inputs, and the script calculates and visualizes the circle centers within the rectangle.

## Features
- Compute centers for circles packed inside a rectangle.
- Flexible inputs for rectangle diagonal coordinates and circle radius.
- Visualize the rectangle and circles using `matplotlib`.

## Requirements
- Python 3.x
- `numpy` library
- `matplotlib` library

You can install the required libraries by running:
```bash
pip install numpy matplotlib
```

## Function Descriptions
- `calculate_circle_centers(x1, y1, x2, y2, r)`<br>
This function computes the centers of the circles inside the rectangle using hexagonal packing. It returns a list of the centers as (x, y) coordinates.
- `plot_circles_and_rectangle(centers, x1, y1, x2, y2, r)`<br>
This function plots the rectangle and the calculated circle centers. It uses the `matplotlib` library to create a visualization where each circle is drawn inside the rectangle with its center marked.
- `main()`<br>
The main function that prompts the user for input, calls the necessary functions to calculate the circle centers, and plots the output.

## Example
# Input
```mathematica
Enter the first diagonal coordinate of the rectangle (x1 y1): 0 0
Enter the second diagonal coordinate of the rectangle (x2 y2): 20 10
Enter the radius of the circles: 1
```
# Output
![Image showing how circles packed in the rectangle](circle-fitting/Figure_1.png)
