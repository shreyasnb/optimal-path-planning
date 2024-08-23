#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <chrono>
#include <algorithm>

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

int main() {
    string filename = "data.csv";  // Replace with your CSV file path

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

        // Find the convex hull and calculate its area
        vector<Point> hull = convex_hull(vertices);
        double hull_area = polygon_area(hull);

        cout << "The minimum area covering all points (Convex Hull Area) is: " << hull_area << endl;

        cout << "Convex Hull vertices:" << endl;
        for (const auto& vertex : hull) {
            cout << "(" << vertex.x << ", " << vertex.y << ")" << endl;
        }

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
