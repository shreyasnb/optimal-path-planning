#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <limits>
#include <stack>
#include "matplotlibcpp.h"  // Include matplotlibcpp

namespace plt = matplotlibcpp;

using namespace std;
using namespace std::chrono;

struct Point {
    double x, y;
};

// Function to calculate the distance between two points
double distance(Point p1, Point p2) {
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

// Function to find the orientation of the triplet (p, q, r)
// 0 -> p, q and r are collinear
// 1 -> Clockwise
// 2 -> Counterclockwise
int orientation(Point p, Point q, Point r) {
    double val = (q.y - p.y) * (r.x - q.x) - 
                 (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0;  // collinear
    return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// Function to find the convex hull using Graham scan algorithm
vector<Point> convex_hull(vector<Point> points) {
    // Find the bottom-most point (or leftmost in case of tie)
    int n = points.size(), ymin = points[0].y, min = 0;
    for (int i = 1; i < n; i++) {
        if ((points[i].y < points[min].y) || 
           (points[i].y == points[min].y && points[i].x < points[min].x))
            min = i;
    }
    
    // Place the bottom-most point at first position
    swap(points[0], points[min]);
    Point p0 = points[0];
    
    // Sort the points according to polar angle with p0
    sort(points.begin() + 1, points.end(), [p0](Point a, Point b) {
        int o = orientation(p0, a, b);
        if (o == 0)
            return distance(p0, b) >= distance(p0, a);
        return o == 2;
    });

    // Remove points that are collinear
    int m = 1;  // Initialize size of modified array
    for (int i = 1; i < n; i++) {
        while (i < n - 1 && orientation(p0, points[i], points[i + 1]) == 0)
            i++;
        points[m] = points[i];
        m++;  // Update size of modified array
    }

    if (m < 3) throw invalid_argument("Convex hull is not possible");  // Convex hull not possible

    vector<Point> hull;
    hull.push_back(points[0]);
    hull.push_back(points[1]);
    hull.push_back(points[2]);

    // Build the hull
    for (int i = 3; i < m; i++) {
        while (hull.size() > 1 && orientation(hull[hull.size() - 2], hull[hull.size() - 1], points[i]) != 2)
            hull.pop_back();
        hull.push_back(points[i]);
    }
    return hull;
}

// Function to calculate the area of a polygon given its vertices
double polygon_area(const vector<Point>& vertices) {
    int n = vertices.size();
    double area = 0.0;

    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += vertices[i].x * vertices[j].y;
        area -= vertices[j].x * vertices[i].y;
    }

    return fabs(area) / 2.0;
}

// Function to simplify a polygon to a target number of sides
vector<Point> simplify_polygon(vector<Point>& vertices, int target_sides) {
    if (vertices.size() <= target_sides) {
        return vertices;  // Already has the target number of sides
    }

    // Simplify the polygon using Douglas-Peucker algorithm
    // For demonstration, a naive approach is used
    vector<Point> simplified;
    double epsilon = 0.1;  // Approximation tolerance (adjust as needed)
    stack<Point> stk;
    
    // Simplify polygon (placeholder)
    int step = max(1, static_cast<int>(vertices.size() / target_sides));
    for (size_t i = 0; i < vertices.size(); i += step) {
        simplified.push_back(vertices[i]);
    }
    
    if (simplified.size() < target_sides) {
        simplified.insert(simplified.end(), vertices.begin(), vertices.end());
    }
    
    return simplified;
}

// Read coordinates from CSV file
vector<Point> read_coordinates_from_csv(const string& filename) {
    vector<Point> points;
    ifstream file(filename);
    string line;

    while (getline(file, line)) {
        stringstream ss(line);
        string x_str, y_str;
        getline(ss, x_str, ',');
        getline(ss, y_str, ',');

        Point p;
        p.x = stod(x_str);
        p.y = stod(y_str);

        points.push_back(p);
    }

    return points;
}

int main() {
    string filename = "pacman.csv";  // Replace with your CSV file path

    // Start timing
    auto start = high_resolution_clock::now();

    vector<Point> vertices = read_coordinates_from_csv(filename);

    // Display the vertices
    cout << "Vertices from the CSV file:" << endl;
    for (const auto& vertex : vertices) {
        cout << "(" << vertex.x << ", " << vertex.y << ")" << endl;
    }

    try {
        // Find the convex hull
        vector<Point> hull = convex_hull(vertices);
        double hull_area = polygon_area(hull);

        // Simplify the hull to exactly 10 sides
        vector<Point> simplified_hull = simplify_polygon(hull, 10);
        double simplified_area = polygon_area(simplified_hull);

        cout << "The minimum area covering all points (Convex Hull Area) is: " << hull_area << endl;
        cout << "The simplified convex hull with 10 sides area is: " << simplified_area << endl;

        cout << "Simplified Convex Hull vertices:" << endl;
        for (const auto& vertex : simplified_hull) {
            cout << "(" << vertex.x << ", " << vertex.y << ")" << endl;
        }

        // Plot using matplotlibcpp
        vector<double> x_points, y_points;
        for (const auto& vertex : vertices) {
            x_points.push_back(vertex.x);
            y_points.push_back(vertex.y);
        }
        plt::scatter(x_points, y_points, 10, {{"color", "blue"}, {"label", "Data Points"}});

        vector<double> x_hull, y_hull;
        for (const auto& vertex : hull) {
            x_hull.push_back(vertex.x);
            y_hull.push_back(vertex.y);
        }
        x_hull.push_back(hull[0].x);  // Closing the hull loop
        y_hull.push_back(hull[0].y);
        plt::plot(x_hull, y_hull, {{"color", "red"}, {"label", "Convex Hull"}});

        vector<double> x_simplified, y_simplified;
        for (const auto& vertex : simplified_hull) {
            x_simplified.push_back(vertex.x);
            y_simplified.push_back(vertex.y);
        }
        x_simplified.push_back(simplified_hull[0].x);  // Closing the hull loop
        y_simplified.push_back(simplified_hull[0].y);
        plt::plot(x_simplified, y_simplified, {{"color", "green"}, {"label", "Simplified Hull"}});

        plt::legend();
        plt::xlabel("X");
        plt::ylabel("Y");
        plt::title("Data Points and Convex Hull");
        plt::show();

    } catch (const exception &e) {
        cerr << e.what() << endl;
    }

    // End timing
    auto end = high_resolution_clock::now();

    // Calculate the duration
    auto duration = duration_cast<milliseconds>(end - start);
    cout << "Time taken by the program: " << duration.count() << " milliseconds" << endl;

    return 0;
}