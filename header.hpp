#ifndef HEADER_H
#define HEADER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <limits>
#include <cmath>

struct Point
{
    double x;
    double y;
};

typedef std::vector<Point> vecP;
typedef std::vector<std::vector<double>> vec2Df;

void findEuclideanMotion(vecP&, vecP&);

#endif  /* HEADER_H */