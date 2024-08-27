#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

struct Point {
    double x, y;
};

// Function to read points from a CSV file
std::vector<Point> readPointsFromCSV(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        return points;
    }

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string x_str, y_str;
        if (std::getline(ss, x_str, ',') && std::getline(ss, y_str, ',')) {
            Point p;
            p.x = std::stod(x_str);
            p.y = std::stod(y_str);
            points.push_back(p);
        }
    }

    file.close();
    return points;
}

// Function to calculate the area of a polygon given its vertices
double polygonArea(const std::vector<Point>& vertices) {
    double area = 0.0;
    int n = vertices.size();
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += vertices[i].x * vertices[j].y;
        area -= vertices[j].x * vertices[i].y;
    }
    return fabs(area) / 2.0;
}

// Function to compute the convex hull using the Graham scan algorithm
std::vector<Point> convexHull(std::vector<Point>& points) {
    int n = points.size();
    std::sort(points.begin(), points.end(), [](Point a, Point b) {
        return a.y < b.y || (a.y == b.y && a.x < b.x);
    });

    std::vector<Point> hull;
    for (int phase = 0; ++phase < 2;) {
        int start = hull.size();
        for (const Point& p : points) {
            while (hull.size() >= start + 2 &&
                   (hull[hull.size() - 1].x - hull[hull.size() - 2].x) * (p.y - hull[hull.size() - 2].y) -
                   (hull[hull.size() - 1].y - hull[hull.size() - 2].y) * (p.x - hull[hull.size() - 2].x) <= 0) {
                hull.pop_back();
            }
            hull.push_back(p);
        }
        hull.pop_back();
        std::reverse(points.begin(), points.end());
    }

    return hull;
}

// Function to simplify a convex hull into a polygon with no more than 10 sides
std::vector<Point> simplifyPolygon(std::vector<Point>& hull, int maxSides) {
    while (hull.size() > maxSides) {
        double minAreaIncrease = std::numeric_limits<double>::max();
        int removeIndex = 0;

        for (int i = 0; i < hull.size(); i++) {
            std::vector<Point> tempHull = hull;
            tempHull.erase(tempHull.begin() + i);
            double area = polygonArea(tempHull);
            if (area < minAreaIncrease) {
                minAreaIncrease = area;
                removeIndex = i;
            }
        }

        hull.erase(hull.begin() + removeIndex);
    }

    return hull;
}

// Function to plot points
void plotPoints(const std::vector<Point>& points, const std::string& label) {
    std::vector<double> x, y;
    for (const auto& p : points) {
        x.push_back(p.x);
        y.push_back(p.y);
    }
    plt::scatter(x, y, 10, {{"label", label}});
}

// Function to plot polygon
void plotPolygon(const std::vector<Point>& polygon, const std::string& label, const std::string& color) {
    std::vector<double> x, y;
    for (const auto& p : polygon) {
        x.push_back(p.x);
        y.push_back(p.y);
    }
    x.push_back(polygon[0].x);  // Closing the polygon
    y.push_back(polygon[0].y);
    plt::plot(x, y, {{"label", label}, {"color", color}});
}

int main() {
    // Read the data points from the CSV
    std::vector<Point> points = readPointsFromCSV("pacman.csv");

    // Step 1: Compute the Convex Hull
    std::vector<Point> hull = convexHull(points);

    // Step 2: Simplify the hull to have at most 10 sides
    std::vector<Point> optimizedPolygon = simplifyPolygon(hull, 10);

    // Step 3: Calculate the area of the optimized polygon
    double optimizedArea = polygonArea(optimizedPolygon);

    // Output the results
    std::cout << "Optimized Polygon Area: " << optimizedArea << std::endl;

    // Plot the points
    plotPoints(points, "Scatter Points");

    // Plot the Convex Hull
    plotPolygon(hull, "Convex Hull", "blue");

    // Plot the Optimized Polygon
    plotPolygon(optimizedPolygon, "Optimized Polygon", "red");

    // Show the plot with legend
    plt::legend();
    plt::show();

    return 0;
}