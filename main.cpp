#include "header.hpp"

void readPointFromFile(vecP& points, std::ifstream& file)
{
    double x = 0;
    double y = 0;
    std::string line;

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        if (!(iss >> x >> y))
            throw std::runtime_error("File read error!");
        if (iss.peek() != EOF)
            throw std::runtime_error("File format error!");
        if (x > std::numeric_limits<double>::max() || x < std::numeric_limits<double>::min() ||
            y > std::numeric_limits<double>::max() || y < std::numeric_limits<double>::min())
            throw std::runtime_error("Point value is not in the double range!");
        points.push_back({x, y});
    }
}

int main() try
{
    std::ifstream fileSetA("input_set_a.txt");
    std::ifstream fileSetB("input_set_b.txt");
    if (!fileSetA.is_open() || !fileSetB.is_open())
        throw std::runtime_error("File open error!");

    // Read setA points from file
    vecP pointsSetA;
    readPointFromFile(pointsSetA, fileSetA);
    fileSetA.close();

    // Read setB points from file
    vecP pointsSetB;
    readPointFromFile(pointsSetB, fileSetB);
    fileSetB.close();

    if (pointsSetA.size() != pointsSetA.size())
        throw std::runtime_error("Sets has different points count");

    findEuclideanMotion(pointsSetA, pointsSetB);
    return 0;
}
catch(const std::exception& e)
{
    std::cerr << e.what() << std::endl;
}