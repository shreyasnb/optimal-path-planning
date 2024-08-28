#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>
#include "gnuplot-iostream.h"

// Function to calculate the Convex Hull using the Gift Wrapping algorithm
std::vector<Eigen::Vector2d> convexHull(const std::vector<Eigen::Vector2d>& points) {
    int n = points.size();
    if (n < 3) return points;

    std::vector<Eigen::Vector2d> hull;

    int leftmost = 0;
    for (int i = 1; i < n; i++)
        if (points[i][0] < points[leftmost][0])
            leftmost = i;

    int p = leftmost, q;
    do {
        hull.push_back(points[p]);
        q = (p + 1) % n;

        for (int i = 0; i < n; i++) {
            if ((points[i] - points[p]).cross(points[q] - points[p]) > 0)
                q = i;
        }

        p = q;

    } while (p != leftmost);

    return hull;
}

// Function to reduce/simplify the convex hull to a 10-sided polygon
std::vector<Eigen::Vector2d> reduceToPolygon(const std::vector<Eigen::Vector2d>& hull_points, int num_sides) {
    int current_sides = hull_points.size();
    std::vector<Eigen::Vector2d> reduced_points;

    if (current_sides == num_sides) {
        return hull_points;
    }

    if (current_sides > num_sides) {
        int step = current_sides / num_sides;
        for (int i = 0; i < num_sides; ++i) {
            reduced_points.push_back(hull_points[(i * step) % current_sides]);
        }
        return reduced_points;
    }

    // Interpolation if fewer points are present
    return interpolatePoints(hull_points, num_sides);
}

// Function to interpolate additional points
std::vector<Eigen::Vector2d> interpolatePoints(const std::vector<Eigen::Vector2d>& hull_points, int num_sides) {
    int current_sides = hull_points.size();
    std::vector<Eigen::Vector2d> new_points;

    for (int i = 0; i < current_sides; ++i) {
        new_points.push_back(hull_points[i]);
        Eigen::Vector2d next_point = hull_points[(i + 1) % current_sides];
        int num_interpolations = (num_sides - current_sides) / current_sides;

        for (int j = 1; j <= num_interpolations; ++j) {
            Eigen::Vector2d interpolated_point = hull_points[i] + (next_point - hull_points[i]) * (j / static_cast<double>(num_interpolations + 1));
            new_points.push_back(interpolated_point);
        }

        if (new_points.size() >= num_sides) {
            break;
        }
    }

    return std::vector<Eigen::Vector2d>(new_points.begin(), new_points.begin() + num_sides);
}

// Check if all points are inside the polygon
bool allInsidePolygon(const std::vector<Eigen::Vector2d>& polygon, const std::vector<Eigen::Vector2d>& points) {
    for (const auto& point : points) {
        bool inside = false;
        int n = polygon.size();
        for (int i = 0, j = n - 1; i < n; j = i++) {
            if (((polygon[i][1] > point[1]) != (polygon[j][1] > point[1])) &&
                (point[0] < (polygon[j][0] - polygon[i][0]) * (point[1] - polygon[i][1]) / (polygon[j][1] - polygon[i][1]) + polygon[i][0])) {
                inside = !inside;
            }
        }
        if (!inside) return false;
    }
    return true;
}

// Function to scale the polygon outwards
std::vector<Eigen::Vector2d> scalePolygon(const std::vector<Eigen::Vector2d>& polygon_points, const std::vector<Eigen::Vector2d>& points, double scale_factor = 1.1) {
    Eigen::Vector2d centroid = Eigen::Vector2d::Zero();
    for (const auto& point : polygon_points) {
        centroid += point;
    }
    centroid /= polygon_points.size();

    std::vector<Eigen::Vector2d> scaled_polygon;
    for (const auto& point : polygon_points) {
        scaled_polygon.push_back((point - centroid) * scale_factor + centroid);
    }

    while (!allInsidePolygon(scaled_polygon, points)) {
        scale_factor += 0.05;
        scaled_polygon.clear();
        for (const auto& point : polygon_points) {
            scaled_polygon.push_back((point - centroid) * scale_factor + centroid);
        }
    }

    return scaled_polygon;
}

int main() {
    std::vector<Eigen::Vector2d> points;

    // Reading the CSV file
    std::ifstream file("pacman.csv");
    std::string line;
    while (std::getline(file, line)) {
        double x, y;
        sscanf(line.c_str(), "%lf,%lf", &x, &y);
        points.emplace_back(x, y);
    }
    file.close();

    // Calculate the convex hull
    std::vector<Eigen::Vector2d> hull_points = convexHull(points);

    // Reduce or expand to form a 10-sided polygon
    std::vector<Eigen::Vector2d> polygon_points = reduceToPolygon(hull_points, 10);

    // Scale the polygon to ensure all points are inside
    std::vector<Eigen::Vector2d> scaled_polygon_points = scalePolygon(polygon_points, points);

    // Plotting using gnuplot
    Gnuplot gp;
    gp << "set title 'Scaled 10-sided Polygon Enclosing All Points'\n";
    gp << "plot '-' with points title 'Data Points', '-' with lines title 'Scaled 10-sided Polygon'\n";
    gp.send1d(points);
    gp.send1d(scaled_polygon_points);

    return 0;
}
