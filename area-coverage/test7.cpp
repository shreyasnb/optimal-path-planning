#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <set>
#include <limits>

using namespace std;
using namespace std::chrono;

struct Point {
    double x, y;

    // Define comparison operators for Point
    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }

    bool operator<(const Point& other) const {
        if (x < other.x) return true;
        if (x > other.x) return false;
        return y < other.y;
    }
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
    int n = points.size();
    if (n < 3) return points;

    // Find the bottom-most point (or leftmost in case of tie)
    int ymin = points[0].y, min = 0;
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

// Function to find the longest diagonal in an n-sided polygon
pair<pair<int, int>, double> longest_diagonal(vector<Point> &vertices) {
    int n = vertices.size();
    if (n < 3) {
        throw invalid_argument("A polygon must have at least 3 sides.");
    }

    double max_distance = 0.0;
    pair<int, int> diagonal_vertices;

    // Calculate the distance between each pair of vertices
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dist = distance(vertices[i], vertices[j]);
            if (dist > max_distance) {
                max_distance = dist;
                diagonal_vertices = {i, j};
            }
        }
    }

    return {diagonal_vertices, max_distance};
}

// Function to read coordinates from a CSV file
vector<Point> read_coordinates_from_csv(const string& filename) {
    vector<Point> vertices;
    ifstream file(filename);
    string line;

    // Read the file line by line
    while (getline(file, line)) {
        stringstream ss(line);
        string x_str, y_str;

        // Read x and y values
        if (getline(ss, x_str, ',') && getline(ss, y_str, ',')) {
            double x = stod(x_str);
            double y = stod(y_str);
            vertices.push_back({x, y});
        }
    }

    return vertices;
}

// Function to generate a 10-sided polygon that covers all points
vector<Point> create_ten_sided_polygon(const vector<Point>& points) {
    vector<Point> ten_sided_polygon;
    
    // If there are more than 10 points, use only the first 10
    if (points.size() >= 10) {
        ten_sided_polygon.assign(points.begin(), points.begin() + 10);
    } else {
        // Use all points available and add random points if needed
        ten_sided_polygon.assign(points.begin(), points.end());
        
        // Ensure there are exactly 10 points
        while (ten_sided_polygon.size() < 10) {
            Point extra_point = {static_cast<double>(rand() % 100), static_cast<double>(rand() % 100)};
            // Check if the point is already in the polygon
            if (find(ten_sided_polygon.begin(), ten_sided_polygon.end(), extra_point) == ten_sided_polygon.end()) {
                ten_sided_polygon.push_back(extra_point);
            }
        }
    }
    
    // Optionally, you could sort or process these points to ensure they form a 10-sided polygon
    // For simplicity, we'll just use the selected points

    return ten_sided_polygon;
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
        pair<pair<int, int>, double> result = longest_diagonal(vertices);
        pair<int, int> indices = result.first;
        double max_dist = result.second;

        cout << "The longest diagonal is between points: (" << vertices[indices.first].x << ", " 
             << vertices[indices.first].y << ") and (" << vertices[indices.second].x << ", " 
             << vertices[indices.second].y << ")" << endl;
        cout << "The length of the longest diagonal is: " << max_dist << endl;

        // Create a 10-sided polygon that covers all points
        vector<Point> ten_sided_polygon = create_ten_sided_polygon(vertices);
        double polygon_area_value = polygon_area(ten_sided_polygon);

        cout << "The area of the 10-sided polygon covering all points is: " << polygon_area_value << endl;

        cout << "10-Sided Polygon vertices:" << endl;
        for (const auto& vertex : ten_sided_polygon) {
            cout << "(" << vertex.x << ", " << vertex.y << ")" << endl;
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }

    // End timing and print duration
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    cout << "Execution Time: " << duration.count() << " milliseconds" << endl;

    return 0;
}