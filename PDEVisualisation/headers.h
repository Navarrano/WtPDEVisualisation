#ifndef HEADERS_H
#define HEADERS_H

#include "../PDESolver/PBspVol.h"

struct ChartData
{
	Matrix3D<Point1D> points;
	double xStart, xEnd, yStart, yEnd, zStart, zEnd, step;
	unsigned nbPts;
};

enum Dimension { U, V, W};

const unsigned possibleDensities[6] = { 2, 5, 10, 20, 50, 100 };
const double possibleStepSizes[6] = { 0.5, 0.2, 0.1, 0.05, 0.02, 0.01 };
const unsigned stepsArraySize = 6;

#endif