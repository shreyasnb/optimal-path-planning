#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <chrono>  // Include the chrono library

using namespace std;
using namespace std::chrono;  // For easier usage of chrono

// Function to calculate the distance between two points
double distance(pair<double, double> p1, pair<double, double> p2) {
    return sqrt(pow(p2.first - p1.first, 2) + pow(p2.second - p1.second, 2));
}

// Function to find the longest diagonal in an n-sided polygon
pair<pair<int, int>, double> longest_diagonal(vector<pair<double, double>> &vertices) {
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
vector<pair<double, double>> read_coordinates_from_csv(const string& filename) {
    vector<pair<double, double>> vertices;
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
            vertices.push_back(make_pair(x, y));
        }
    }

    return vertices;
}

int main() {
    string filename = "data.csv";  // Replace with your CSV file path

    // Start timing
    auto start = high_resolution_clock::now();

    vector<pair<double, double>> vertices = read_coordinates_from_csv(filename);

    // Display the vertices
    cout << "Vertices from the CSV file:" << endl;
    for (const auto& vertex : vertices) {
        cout << "(" << vertex.first << ", " << vertex.second << ")" << endl;
    }

    try {
        pair<pair<int, int>, double> result = longest_diagonal(vertices);
        pair<int, int> indices = result.first;
        double max_dist = result.second;

        cout << "The longest diagonal is between points: (" << vertices[indices.first].first << ", " 
             << vertices[indices.first].second << ") and (" << vertices[indices.second].first << ", " 
             << vertices[indices.second].second << ")" << endl;
        cout << "The length of the longest diagonal is: " << max_dist << endl;
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
