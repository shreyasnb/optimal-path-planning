#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <random>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;
using namespace std;
using namespace std::chrono;

struct Point {
    double x, y;
};

// Function to calculate the cross product of two vectors OA and OB
double crossProduct(Point O, Point A, Point B) {
    return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}

// Function to compute the Convex Hull using the Monotone Chain algorithm
vector<Point> convexHull(vector<Point> points) {
    int n = points.size(), k = 0;
    if (n <= 3) return points;

    vector<Point> hull(2 * n);

    // Sort points lexicographically
    sort(points.begin(), points.end(), [](Point a, Point b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });

    // Build the lower hull
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && crossProduct(hull[k - 2], hull[k - 1], points[i]) <= 0) k--;
        hull[k++] = points[i];
    }

    // Build the upper hull
    for (int i = n - 2, t = k + 1; i >= 0; --i) {
        while (k >= t && crossProduct(hull[k - 2], hull[k - 1], points[i]) <= 0) k--;
        hull[k++] = points[i];
    }

    hull.resize(k - 1);
    return hull;
}

// Function to compute the area of a polygon given its vertices
double polygonArea(const vector<Point> &polygon) {
    double area = 0.0;
    int n = polygon.size();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += polygon[i].x * polygon[j].y;
        area -= polygon[i].y * polygon[j].x;
    }
    return abs(area) / 2.0;
}

// Function to simplify a polygon to a maximum of 10 sides
vector<Point> simplifyPolygon(const vector<Point> &polygon, int maxSides) {
    vector<Point> simplified = polygon;

    while (simplified.size() > maxSides) {
        double minLoss = numeric_limits<double>::max();
        int indexToRemove = -1;

        for (int i = 0; i < simplified.size(); ++i) {
            int prev = (i == 0) ? simplified.size() - 1 : i - 1;
            int next = (i + 1) % simplified.size();

            double loss = abs(crossProduct(simplified[prev], simplified[i], simplified[next]));

            if (loss < minLoss) {
                minLoss = loss;
                indexToRemove = i;
            }
        }

        simplified.erase(simplified.begin() + indexToRemove);
    }

    return simplified;
}

// Function to calculate the distance between two points
double distance(Point p1, Point p2) {
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

// Function to find the longest diagonal in a set of points
pair<pair<int, int>, double> longest_diagonal(vector<Point> &points) {
    int n = points.size();
    if (n < 3) {
        throw invalid_argument("A polygon must have at least 3 sides.");
    }

    double max_distance = 0.0;
    pair<int, int> diagonal_vertices;

    // Calculate the distance between each pair of points
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dist = distance(points[i], points[j]);
            if (dist > max_distance) {
                max_distance = dist;
                diagonal_vertices = {i, j};
            }
        }
    }

    return {diagonal_vertices, max_distance};
}

// Function to load points from a CSV file
vector<Point> loadPointsFromCSV(const string &filename) {
    vector<Point> points;
    ifstream file(filename);
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            double x, y;
            sscanf(line.c_str(), "%lf,%lf", &x, &y);
            points.push_back({x, y});
        }
        file.close();
    }
    return points;
}

// Function to plot the points and polygons using Matplotlib for C++
void plotPointsAndPolygons(const vector<Point> &points, const vector<Point> &hull, const vector<Point> &polygon) {
    // Plot points
    for (const auto &p : points) {
        std::vector<double> x = {p.x};
        std::vector<double> y = {p.y};
        plt::scatter(x, y, 10.0);
    }

    // Plot Convex Hull
    vector<double> hx, hy;
    for (const auto &p : hull) {
        hx.push_back(p.x);
        hy.push_back(p.y);
    }
    hx.push_back(hx[0]); // Close the hull
    hy.push_back(hy[0]);
    plt::plot(hx, hy, "r-");

    // Plot 10-sided Polygon
    vector<double> px, py;
    for (const auto &p : polygon) {
        px.push_back(p.x);
        py.push_back(p.y);
    }
    px.push_back(px[0]); // Close the polygon
    py.push_back(py[0]);
    plt::plot(px, py, "g-");

    // Display the plot
    plt::show();
}

int main() {
    // Start timing
    auto start = high_resolution_clock::now();

    // Load points from CSV
    vector<Point> points = loadPointsFromCSV("data2.csv");
    cout << "Number of points read from CSV: " << points.size() << endl;

    // Compute the Convex Hull
    vector<Point> hull = convexHull(points);

    // Simplify the Convex Hull to a 10-sided polygon
    vector<Point> bestPolygon = simplifyPolygon(hull, 10);

    // Calculate the area of the Convex Hull and the 10-sided polygon
    double hullArea = polygonArea(hull);
    double polygonAreaVal = polygonArea(bestPolygon);

    // Find the longest diagonal in the original points
    pair<pair<int, int>, double> diagonal = longest_diagonal(points);
    cout << "The longest diagonal is between points: (" << points[diagonal.first.first].x << ", " 
         << points[diagonal.first.first].y << ") and (" << points[diagonal.first.second].x << ", " 
         << points[diagonal.first.second].y << ") with a distance of " << diagonal.second << endl;

    // Output the results
    cout << "Convex Hull Area: " << hullArea << endl;
    cout << "10-sided Polygon Area: " << polygonAreaVal << endl;

    // Stop timing
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Time taken to run the program: " << duration.count() << " milliseconds" << endl;

    // Plot the points, Convex Hull, and 10-sided polygon
    plotPointsAndPolygons(points, hull, bestPolygon);

    return 0;
}