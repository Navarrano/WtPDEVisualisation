#ifndef HEADERS_H
#define HEADERS_H

#include "../PDESolver/PBspVol.h"
#include <sstream>
#include <iomanip>



struct ChartData
{
	Matrix3D<Point1D> points;
	double xStart, xEnd, yStart, yEnd, zStart, zEnd, step;
	unsigned xSize, ySize, zSize;
	unsigned totalNbPts;
};

struct Limits
{
	Limits() {};
	Limits(ChartData& data) : xLeftLimit(data.xStart), xRightLimit(data.xEnd), yLeftLimit(data.xStart), yRightLimit(data.yEnd), zLeftLimit(data.zStart), zRightLimit(data.zEnd) {};
	Limits(double xL, double xR, double yL, double yR, double zL, double zR) : xLeftLimit(xL), xRightLimit(xR), yLeftLimit(yL), yRightLimit(yR), zLeftLimit(zL), zRightLimit(zR) {};
	double xLeftLimit, xRightLimit;
	double yLeftLimit, yRightLimit;
	double zLeftLimit, zRightLimit;
};

struct SurfaceOrders
{
	SurfaceOrders() {};
	SurfaceOrders(short x, short y, short z) : xOrder(x), yOrder(y), zOrder(z) {};
	short xOrder, yOrder, zOrder;
};

enum Dimension { U, V, W};

enum DisplayRate { Every = 1, EverySecond = 2, EveryFifth = 5, EveryTenth = 10};

const unsigned possibleDensities[6] = { 2, 5, 10, 20, 50, 100 };
const double possibleStepSizes[6] = { 0.5, 0.2, 0.1, 0.05, 0.02, 0.01 };
const unsigned stepsArraySize = 6;

#endif