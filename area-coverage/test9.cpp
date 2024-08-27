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

    Point(double x, double y) : x(x), y(y) {}

    bool operator<(const Point& other) const {
        if (x == other.x) return y < other.y;
        return x < other.x;
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

// Function to generate new points to create a 10-sided polygon
vector<Point> add_points_to_polygon(vector<Point>& vertices, int required_sides) {
    vector<Point> new_polygon;
    int n = vertices.size();

    if (n >= required_sides) {
        // If there are enough vertices, just return the polygon
        new_polygon = vertices;
    } else {
        // Generate new points to make up the difference
        set<Point> unique_points(vertices.begin(), vertices.end());
        new_polygon = vertices;
        
        int points_needed = required_sides - n;
        int steps = n;
        
        while (points_needed > 0) {
            for (int i = 0; i < n; ++i) {
                if (points_needed <= 0) break;
                Point current = vertices[i];
                Point next = vertices[(i + 1) % n];
                Point mid = {(current.x + next.x) / 2, (current.y + next.y) / 2};
                if (unique_points.find(mid) == unique_points.end()) {
                    new_polygon.push_back(mid);
                    unique_points.insert(mid);
                    points_needed--;
                }
            }
        }
    }
    
    // Ensure the polygon is closed
    if (new_polygon.size() > 1 && new_polygon.front() != new_polygon.back()) {
        new_polygon.push_back(new_polygon.front());
    }

    return new_polygon;
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

// Updated function to create a 10-sided polygon
vector<Point> create_ten_sided_polygon(const vector<Point>& points) {
    // Start with the convex hull
    vector<Point> hull = convex_hull(points);

    // Add new points to form a 10-sided polygon
    vector<Point> ten_sided_polygon = add_points_to_polygon(hull, 10);

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
        // Create a 10-sided polygon that covers all points
        vector<Point> ten_sided_polygon = create_ten_sided_polygon(vertices);
        double polygon_area_value = polygon_area(ten_sided_polygon);

        cout << "The area of the 10-sided polygon covering all points is: " << polygon_area_value << endl;

        cout << "10-Sided Polygon vertices:" << endl;
        for (const auto& vertex : ten_sided_polygon) {
            cout << "(" << vertex.x << ", " << vertex.y << ")" << endl;
        }

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
