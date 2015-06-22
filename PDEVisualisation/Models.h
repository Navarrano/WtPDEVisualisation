#ifndef MODELS_H
#define MODELS_H

#include "../PDESolver/PBspSurf.h"
#include "../PDESolver/PBspVol.h"
#include "headers.h"

#include <Wt/WAbstractTableModel>
#include <Wt/WStandardItemModel>
#include <Wt/Chart/WStandardColorMap>

using namespace Wt;

class SurfaceData : public WStandardItemModel
{
public:
	SurfaceData(Matrix<Point1D> points, double xStart, double yStart, double step, WObject *parent = 0);
	boost::any data(const WModelIndex& index, int role = DisplayRole) const;
private:
	const Matrix<Point1D> points_;
	const double xStart_, yStart_, step_;
};

class VolumeData : public Wt::WStandardItemModel {
public:
	//VolumeData(Matrix3D<Point1D> points, double xStart, double yStart, double zStart, double step, Limits limits, Wt::WObject *parent = 0);
	VolumeData(ChartData& chartData, WObject *parent = 0);
	virtual ~VolumeData();
	boost::any data(const Wt::WModelIndex& index,
		int role = Wt::DisplayRole) const;
	double getMin() { return min_; };
	double getMax() { return max_; };
	void setXRate(DisplayRate rate) { xRate_ = rate; };
	void setYRate(DisplayRate rate) { yRate_ = rate; };
	void setZRate(DisplayRate rate) { zRate_ = rate; };
	void setPointSize(short size) { pointSize_ = size; };
	void setSkippedPointSize(short size) { skippedPointSize_ = size; };
	void setLimits(Limits limits) { limits_ = limits; };
private:
	void findMinAndMax();
	bool checkXLimits(double x) const;
	bool checkYLimits(double y) const;
	bool checkZLimits(double z) const;

	const unsigned xSize_, ySize_, zSize_;
	const double xStart_, yStart_, zStart_, step_;
	double min_, max_;
	short pointSize_, skippedPointSize_;
	DisplayRate xRate_, yRate_, zRate_;
	Matrix3D<Point1D> points_;
	Chart::WStandardColorMap *colorMap_;
	Limits limits_;
};

#endif
