#ifndef MODELS_H
#define MODELS_H

#include "../PDESolver/PBspSurf.h"
#include "../PDESolver/PBspVol.h"

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

class SurfaceData2 : public WStandardItemModel
{
public:
	SurfaceData2(FBspSurf& surf, unsigned nbXpts, unsigned nbYpts, WObject *parent = 0);

	boost::any data(const WModelIndex& index, int role = DisplayRole) const;
private:
	Matrix<double> data_;
	const double xStart_, xEnd_, yStart_, yEnd_, xStep_, yStep_;
};

class VolumeData : public Wt::WStandardItemModel {
public:
	VolumeData(Matrix3D<Point1D> points, double xStart, double yStart, double zStart, double step, Wt::WObject *parent = 0);
	virtual ~VolumeData();
	boost::any data(const Wt::WModelIndex& index,
		int role = Wt::DisplayRole) const;
	double getMin() { return min_; };
	double getMax() { return max_; };
private:
	void findMinAndMax();

	const unsigned nbPts_;
	const double xStart_, yStart_, zStart_, step_;
	double min_, max_;
	Matrix3D<Point1D> points_;
	Chart::WStandardColorMap *colorMap_;
};

#endif
