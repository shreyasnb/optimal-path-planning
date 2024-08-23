#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <random>

using namespace std;

// Function to generate a random pair of doubles within the specified range
std::pair<double, double> generate_random_pair(double min, double max) {
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator

    // Define distributions for each range
    std::uniform_real_distribution<> dis(min, max);

    // Generate the random pair
    return std::make_pair(dis(gen), dis(gen));
}

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

int main() {
    int n = 0;
    cout << "Enter the number of sides: ";
    cin >> n;

    vector<pair<double, double>> vertices;

    // Generate random vertices
    for (int i = 0; i < n; ++i) {
        vertices.push_back(generate_random_pair(0, 100));
    }

    // Display the generated vertices
    cout << "Generated vertices:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "(" << vertices[i].first << ", " << vertices[i].second << ")" << endl;
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

    return 0;
}
