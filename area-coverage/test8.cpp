#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <chrono>
#include <limits>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;
using namespace std::chrono;

struct Point {

    double x, y;

    // Define the equality operator
    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }


};

// Function to read coordinates from a CSV file
vector<Point> read_coordinates_from_csv(const string& filename) {
    vector<Point> points;
    ifstream file(filename);
    string line;

    while (getline(file, line)) {
        stringstream ss(line);
        string item;
        Point p;

        getline(ss, item, ',');
        p.x = stod(item);
        getline(ss, item, ',');
        p.y = stod(item);

        points.push_back(p);
    }

    return points;
}

// Function to calculate the distance between two points
double distance(const Point& p1, const Point& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Helper function to calculate the distance from a point to a line segment
double distance_to_line_segment(const Point& p, const Point& v, const Point& w) {
    double l2 = distance(v, w);
    if (l2 == 0.0) return distance(p, v);
    double t = max(0.0, min(1.0, ((p.x - v.x) * (w.x - v.x) + (p.y - v.y) * (w.y - v.y)) / l2));
    Point projection = {v.x + t * (w.x - v.x), v.y + t * (w.y - v.y)};
    return distance(p, projection);
}

// Function to calculate the cross product of vectors OA and OB
double cross_product(const Point& O, const Point& A, const Point& B) {
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// Function to calculate the convex hull of a set of points using Andrew's monotone chain algorithm
vector<Point> convex_hull(vector<Point> points) {
    if (points.size() < 3) throw runtime_error("At least 3 points are required to form a convex hull.");

    // Sort the points lexicographically
    sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });

    // Build the lower hull
    vector<Point> lower;
    for (const auto& p : points) {
        while (lower.size() >= 2 && cross_product(lower[lower.size() - 2], lower.back(), p) <= 0) {
            lower.pop_back();
        }
        lower.push_back(p);
    }

    // Build the upper hull
    vector<Point> upper;
    for (auto it = points.rbegin(); it != points.rend(); ++it) {
        const auto& p = *it;
        while (upper.size() >= 2 && cross_product(upper[upper.size() - 2], upper.back(), p) <= 0) {
            upper.pop_back();
        }
        upper.push_back(p);
    }

    // Remove the last point of each half because it is repeated at the beginning of the other half
    lower.pop_back();
    upper.pop_back();

    // Concatenate lower and upper hulls to get the convex hull
    lower.insert(lower.end(), upper.begin(), upper.end());

    return lower;
}

// Function to create a 10-sided polygon that includes all points
vector<Point> create_ten_sided_polygon(const vector<Point>& points) {
    // First, compute the convex hull
    vector<Point> hull = convex_hull(points);

    // If the hull already has 10 vertices, return it
    if (hull.size() == 10) {
        return hull;
    }
    // If the hull has more than 10 vertices, downsample to 10
    else if (hull.size() > 10) {
        // Downsample using a simple method (e.g., removing every nth point)
        vector<Point> ten_sided_polygon;
        int n = hull.size() / 10;  // Simple downsampling step
        for (size_t i = 0; i < hull.size(); i += n) {
            ten_sided_polygon.push_back(hull[i]);
            if (ten_sided_polygon.size() == 10) break;
        }
        return ten_sided_polygon;
    }
    // If the hull has fewer than 10 vertices, add additional points
    else {
        vector<Point> ten_sided_polygon = hull;

        while (ten_sided_polygon.size() < 10) {
            // Find the point furthest from the current polygon and add it
            double max_distance = -numeric_limits<double>::infinity();
            size_t furthest_point_index = 0;

            for (size_t i = 0; i < points.size(); ++i) {
                // Check if the point is already in the polygon
                if (find(ten_sided_polygon.begin(), ten_sided_polygon.end(), points[i]) != ten_sided_polygon.end()) {
                    continue;
                }

                // Find the minimum distance from this point to any edge of the polygon
                double min_dist_to_polygon = numeric_limits<double>::infinity();
                for (size_t j = 0; j < ten_sided_polygon.size(); ++j) {
                    size_t next_j = (j + 1) % ten_sided_polygon.size();
                    double dist = distance_to_line_segment(points[i], ten_sided_polygon[j], ten_sided_polygon[next_j]);
                    min_dist_to_polygon = min(min_dist_to_polygon, dist);
                }

                // Update the furthest point
                if (min_dist_to_polygon > max_distance) {
                    max_distance = min_dist_to_polygon;
                    furthest_point_index = i;
                }
            }

            // Add the furthest point to the polygon
            ten_sided_polygon.push_back(points[furthest_point_index]);
        }

        return ten_sided_polygon;
    }
}

// Function to calculate the area of a polygon
double polygon_area(const vector<Point>& poly) {
    double area = 0;
    int n = poly.size();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += (poly[i].x * poly[j].y) - (poly[j].x * poly[i].y);
    }
    return fabs(area) / 2.0;
}

// Function to plot points, convex hull, and 10-sided polygon
void plot(const vector<Point>& points, const vector<Point>& hull, const vector<Point>& ten_sided_polygon) {
    vector<double> x_points, y_points;
    vector<double> x_hull, y_hull;
    vector<double> x_polygon, y_polygon;

    for (const auto& p : points) {
        x_points.push_back(p.x);
        y_points.push_back(p.y);
    }

    for (const auto& p : hull) {
        x_hull.push_back(p.x);
        y_hull.push_back(p.y);
    }

    for (const auto& p : ten_sided_polygon) {
        x_polygon.push_back(p.x);
        y_polygon.push_back(p.y);
    }

    // Close the hull and polygon by connecting the last point to the first
    if (!x_hull.empty()) {
        x_hull.push_back(x_hull.front());
        y_hull.push_back(y_hull.front());
    }

    if (!x_polygon.empty()) {
        x_polygon.push_back(x_polygon.front());
        y_polygon.push_back(y_polygon.front());
    }

    // Plot the dataset
    plt::scatter(x_points, y_points, 10.0, {{"label", "Data Points"}});

    // Plot the convex hull
    plt::plot(x_hull, y_hull, {{"label", "Convex Hull"}, {"color", "r"}});

    // Plot the 10-sided polygon
    plt::plot(x_polygon, y_polygon, {{"label", "10-Sided Polygon"}, {"color", "g"}});

    // Add legend and show plot
    plt::legend();
    plt::show();
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
        // Create a 10-sided polygon that covers all points
        vector<Point> ten_sided_polygon = create_ten_sided_polygon(vertices);
        double polygon_area_value = polygon_area(ten_sided_polygon);

        cout << "The area of the 10-sided polygon covering all points is: " << polygon_area_value << endl;

        cout << "10-Sided Polygon vertices:" << endl;
        for (const auto& vertex : ten_sided_polygon) {
            cout << "(" << vertex.x << ", " << vertex.y << ")" << endl;
        }

        // Create the convex hull
        vector<Point> hull = convex_hull(vertices);

        // Plot the dataset, convex hull, and 10-sided polygon
        plot(vertices, hull, ten_sided_polygon);

    } catch (const exception &e) {
        cerr << e.what() << endl;
    }

    // End timing
    auto end = high_resolution_clock::now();

    // Calculate the duration
    auto duration = duration_cast<milliseconds>(end - start);

    cout << "Execution Time: " << duration.count() << " milliseconds" << endl;

    return 0;
}
