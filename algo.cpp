#include "header.hpp"

void computeMean(Point& mean, const vecP& points)
{
    for (const auto& point : points)
    {
        mean.x += point.x;
        mean.y += point.y;
    }

    mean.x /= points.size();
    mean.y /= points.size();
}

void centerData(const Point& mean, vecP& points)
{
    for (auto& point : points)
    {
        point.x -= mean.x;
        point.y -= mean.y;
    }
}

void computeCovMatrix(vec2Df& covMatrix, const vecP& points, const Point& mean)
{
    for (const auto& point : points)
    {
        covMatrix[0][0] += (point.x - mean.x) * (point.x - mean.x);
        covMatrix[0][1] += (point.x - mean.x) * (point.y - mean.y);
        covMatrix[1][0] += (point.y - mean.y) * (point.x - mean.x);
        covMatrix[1][1] += (point.y - mean.y) * (point.y - mean.y);
    }

    double scale = 1.0 / (points.size() - 1);
    for (auto& row : covMatrix)
    {
        for (auto& value : row)
            value *= scale;
    }
}

void computeEigen(const vec2Df& covMatrix, std::vector<double>& eigenvalues,
                  vec2Df& eigenvectors)
{
    double a = covMatrix[0][0];
    double b = covMatrix[0][1];
    double c = covMatrix[1][0];
    double d = covMatrix[1][1];

    // Calculate eigen values
    double discriminant = std::sqrt((a + d) * (a + d) - 4 * (a * d - b * c));
	if (discriminant < 0)
		throw std::runtime_error("Dicriminant for eigen values is less than 0!");
    double lambda1 = (a + d + discriminant) / 2.0;
    double lambda2 = (a + d - discriminant) / 2.0;
    eigenvalues = {lambda1, lambda2};

    // Calculate eigen vectors
    if (b != 0)
	{
        double v1x = lambda1 - d;
        double v1y = b;
        double norm1 = std::sqrt(v1x * v1x + v1y * v1y);
        v1x /= norm1;
        v1y /= norm1;
        eigenvectors.push_back({v1x, v1y});

        double v2x = lambda2 - d;
        double v2y = b;
        double norm2 = std::sqrt(v2x * v2x + v2y * v2y);
        v2x /= norm2;
        v2y /= norm2;
        eigenvectors.push_back({v2x, v2y});
    }
	else if (c != 0)
	{
        double v1x = c;
        double v1y = lambda1 - a;
        double norm1 = std::sqrt(v1x * v1x + v1y * v1y);
        v1x /= norm1;
        v1y /= norm1;
        eigenvectors.push_back({v1x, v1y});

        double v2x = c;
        double v2y = lambda2 - a;
        double norm2 = std::sqrt(v2x * v2x + v2y * v2y);
        v2x /= norm2;
        v2y /= norm2;
        eigenvectors.push_back({v2x, v2y});
    }
	else // b and c is 0, so matrix not have eigen vectors
        throw std::runtime_error("The matrix does not have eigen vectors.");
}

void applyTransform(const vecP& sourcePoints, vecP& transformedPoints, 
						const std::vector<double>& eigenVector, const double scale)
{
    transformedPoints.clear();

    for (const auto& point : sourcePoints)
	{
        double transformedX = scale * (eigenVector[0] * point.x - eigenVector[1] * point.y);
        double transformedY = scale * (eigenVector[1] * point.x + eigenVector[0] * point.y);
        transformedPoints.push_back({transformedX, transformedY});
    }
}

double evalTransformQuality(const vecP& pointsA, const vecP& transformedPointsB)
{
    double sumSquaredDistances = 0.0;
    double sumAbsoluteDistances = 0.0;

    for (size_t i = 0; i < pointsA.size(); ++i)
    {
        double dx = pointsA[i].x - transformedPointsB[i].x;
        double dy = pointsA[i].y - transformedPointsB[i].y;
        sumSquaredDistances += dx * dx + dy * dy;
        
        double abs_dx = std::abs(dx);
        double abs_dy = std::abs(dy);
        sumAbsoluteDistances += abs_dx + abs_dy;
    }
    
    double rmsd = std::sqrt(sumSquaredDistances / pointsA.size());
    double mae = sumAbsoluteDistances / pointsA.size();
    double overallQuality = (rmsd + mae) / 2.0;
	return overallQuality;
}

void findBestOverallQuality(const vecP& setA, const vecP& setB, std::vector<double>& mainEigenVectorA)
{
	double overallQuality = 0.0;
    double bestQuality = std::numeric_limits<double>::infinity();
	std::vector<Point> transformedSetB;

    for (double scale = 0.0; scale <= 10; scale += 0.05)
	{
		// Apply the transformation to one of the sets of points
    	applyTransform(setB, transformedSetB, mainEigenVectorA, scale);

		// Check the quality of transformation
    	overallQuality = evalTransformQuality(setA, transformedSetB);
        if (overallQuality < bestQuality)
            bestQuality = overallQuality;
    }
    std::cout << "Overall Quality (RMSD + MAE): " << bestQuality << std::endl;
}

void findEuclideanMotion(vecP& setA, vecP& setB)
{
    Point meanSetA = {0.0, 0.0};
    Point meanSetB = {0.0, 0.0};

    // Calculating the average value of sets
    computeMean(meanSetA, setA);
    computeMean(meanSetB, setB);

    // Center data of sets
    centerData(meanSetA, setA);
    centerData(meanSetB, setB);

    // Compute covariance matrices
    vec2Df covMatrixA(2, std::vector<double>(2, 0.0));
    vec2Df covMatrixB(2, std::vector<double>(2, 0.0));

    computeCovMatrix(covMatrixA, setA, meanSetA);
    computeCovMatrix(covMatrixB, setB, meanSetB);

	// Calculating eigen values and eigen vectors for covariance matrices
	vec2Df eigenVectorsA;
	vec2Df eigenVectorsB;
	std::vector<double> eigenValuesA;
	std::vector<double> eigenValuesB;

	computeEigen(covMatrixA, eigenValuesA, eigenVectorsA);
	computeEigen(covMatrixB, eigenValuesB, eigenVectorsB);

	// Choosing the main eigen vector
    std::vector<double> mainEigenVectorA = eigenVectorsA[0];

	// Find best transformation
	findBestOverallQuality(setA, setB, mainEigenVectorA);
}