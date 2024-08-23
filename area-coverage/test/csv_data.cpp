#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>  // For std::pair

int main() {
    std::ifstream file("./data.csv");  // Replace with your CSV file path
    std::string line;
    std::vector<std::pair<double, double>> coordinates;

    // Read the file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        double x, y;

        // Read the first value (x coordinate)
        if (std::getline(ss, value, ',')) {
            x = std::stod(value);
        }

        // Read the second value (y coordinate)
        if (std::getline(ss, value, ',')) {
            y = std::stod(value);
        }

        // Store the x and y coordinates as a pair
        coordinates.push_back(std::make_pair(x, y));
    }

    // Process the data (print the coordinates)
    for (const auto& coord : coordinates) {
        std::cout << "Coordinate: (" << coord.first << ", " << coord.second << ")" << std::endl;
    }

    return 0;
}
