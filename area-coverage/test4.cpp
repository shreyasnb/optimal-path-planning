#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <chrono>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;
using namespace std;

struct Point {
    double x, y;
};

// Function to calculate the cross product of vectors OA and OB
double cross_product(const Point& O, const Point& A, const Point& B) {
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// Function to compute the convex hull of a set of points
vector<Point> convex_hull(vector<Point>& points) {
    sort(points.begin(), points.end(), [](Point a, Point b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });

    vector<Point> hull;

    // Build lower hull
    for (const Point& p : points) {
        while (hull.size() >= 2 && cross_product(hull[hull.size() - 2], hull.back(), p) <= 0) {
            hull.pop_back();
        }
        hull.push_back(p);
    }

    // Build upper hull
    for (int i = points.size() - 2, t = hull.size() + 1; i >= 0; --i) {
        while (hull.size() >= t && cross_product(hull[hull.size() - 2], hull.back(), points[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(points[i]);
    }

    hull.pop_back();
    return hull;
}

// Function to compute the area of a polygon
double polygon_area(const vector<Point>& polygon) {
    double area = 0.0;
    int n = polygon.size();

    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += polygon[i].x * polygon[j].y;
        area -= polygon[j].x * polygon[i].y;
    }

    return fabs(area) / 2.0;
}

// Function to reduce the number of vertices of a polygon
vector<Point> reduce_vertices(const vector<Point>& polygon, int max_vertices) {
    if (polygon.size() <= max_vertices) return polygon;

    vector<Point> reduced_polygon = polygon;

    // Simplified approach: reduce by merging consecutive points
    while (reduced_polygon.size() > max_vertices) {
        double min_increase = 1e9;
        int min_index = -1;

        for (int i = 0; i < reduced_polygon.size(); ++i) {
            int prev = (i - 1 + reduced_polygon.size()) % reduced_polygon.size();
            int next = (i + 1) % reduced_polygon.size();

            vector<Point> temp_polygon = reduced_polygon;
            temp_polygon.erase(temp_polygon.begin() + i);

            double area_before = polygon_area(reduced_polygon);
            double area_after = polygon_area(temp_polygon);

            double area_increase = area_after - area_before;

            if (area_increase < min_increase) {
                min_increase = area_increase;
                min_index = i;
            }
        }

        if (min_index != -1) {
            reduced_polygon.erase(reduced_polygon.begin() + min_index);
        }
    }

    return reduced_polygon;
}

// Function to read points from a CSV file
vector<Point> read_points_from_csv(const string& filename) {
    vector<Point> points;
    ifstream file(filename);
    string line;

    while (getline(file, line)) {
        double x, y;
        sscanf(line.c_str(), "%lf,%lf", &x, &y);
        points.push_back({x, y});
    }

    return points;
}

int main() {
    auto start_time = chrono::high_resolution_clock::now();

    // Read the scatter points from data.csv
    vector<Point> points = read_points_from_csv("data.csv");
    cout << "Points read: " << points.size() << endl;

    // Compute the convex hull
    vector<Point> hull = convex_hull(points);
    cout << "Convex hull vertices: " << hull.size() << endl;

    // Reduce vertices if needed
    int max_vertices = 10;
    vector<Point> optimized_polygon = reduce_vertices(hull, max_vertices);
    cout << "Optimized polygon vertices: " << optimized_polygon.size() << endl;

    // Calculate areas
    double hull_area = polygon_area(hull);
    double optimized_area = polygon_area(optimized_polygon);
    cout << "Convex hull area: " << hull_area << endl;
    cout << "Optimized polygon area: " << optimized_area << endl;
    cout << "Excess area: " << optimized_area - hull_area << endl;

    // Plotting
    plt::plot({points.begin(), points.end()}, "ro");
    vector<double> hull_x, hull_y, opt_x, opt_y;

    for (const auto& p : hull) {
        hull_x.push_back(p.x);
        hull_y.push_back(p.y);
    }

    for (const auto& p : optimized_polygon) {
    opt_x.push_back(p.x);
    opt_y.push_back(p.y);
    }

    // Close the polygons by connecting the last point to the first
    hull_x.push_back(hull_x[0]);
    hull_y.push_back(hull_y[0]);
    opt_x.push_back(opt_x[0]);
    opt_y.push_back(opt_y[0]);

    plt::plot(hull_x, hull_y, "b-");  // Convex hull
    plt::plot(opt_x, opt_y, "g-");    // Optimized polygon
    plt::show();

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << "Time taken: " << elapsed.count() << " seconds" << endl;

    return 0;
    }