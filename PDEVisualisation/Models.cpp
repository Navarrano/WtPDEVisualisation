#include "Models.h"

#include "Wt/WColor"

SurfaceData::SurfaceData(Matrix<Point1D> points, double xStart, double yStart, double step, WObject *parent)
: WStandardItemModel(points.GetNumRows() + 1, points.GetNumCols() + 1, parent), points_(points), xStart_(xStart), yStart_(yStart), step_(step)
{}

boost::any SurfaceData::data(const WModelIndex& index, int role) const
{
	if (role != DisplayRole) {
		return WStandardItemModel::data(index, role);
	}
	if (index.row() == 0)
	{
		if (index.column() == 0)
		{
			return 0.0;
		}
		return yStart_ + (index.column() - 1)*step_;
	}
	if (index.column() == 0)
	{
		if (index.row() == 0)
		{
			return 0.0;
		}
		return xStart_ + (index.row() - 1)*step_;
	}
	return points_[index.row() - 1][index.column() - 1].GetX();
}

VolumeData::VolumeData(ChartData& chartData, Wt::WObject *parent) : WStandardItemModel(chartData.totalNbPts, 5, parent), points_(chartData.points),
xStart_(chartData.xStart), yStart_(chartData.yStart), zStart_(chartData.zStart), step_(chartData.step),
min_(INT_MAX), max_(INT_MIN), xSize_(chartData.xSize), ySize_(chartData.ySize), zSize_(chartData.zSize),
xRate_(DisplayRate::Every), yRate_(DisplayRate::Every), zRate_(DisplayRate::Every), pointSize_(6), skippedPointSize_(1), limits_(chartData)
{
	findMinAndMax();
	
	/*std::vector<Chart::WStandardColorMap::Pair> colors;
	colors.push_back(Chart::WStandardColorMap::Pair(min_, WColor(49, 130, 189)));
	colors.push_back(Chart::WStandardColorMap::Pair((min_ + max_) / 2.0, WColor(158, 202, 225)));
	colors.push_back(Chart::WStandardColorMap::Pair(max_, WColor(222, 235, 247)));
	colorMap_ = new Chart::WStandardColorMap(min_, max_, colors, true);*/

	colorMap_ = new Chart::WStandardColorMap(min_, max_, true);
}

void VolumeData::findMinAndMax()
{
	for (int k = 0; k < points_.GetNumPols(); k++)
	{
		for (int i = 0; i < points_.GetNumPols(); i++)
		{
			for (int j = 0; j < points_.GetNumCols(); j++)
			{
				double p = points_[k][i][j].GetX();
				if (p > max_)
				{
					max_ = p;
				}
				else if (p < min_)
				{
					min_ = p;
				}
			}
		}
	}
}

boost::any VolumeData::data(const Wt::WModelIndex& index, int role) const
{
	if (role != DisplayRole) {
		return WStandardItemModel::data(index, role);
	}
	unsigned xIndex = ((index.row() % (xSize_ * ySize_)) / ySize_);
	unsigned yIndex = (index.row() % ySize_);
	unsigned zIndex = (index.row() / (xSize_ * ySize_));
	double x = xStart_ + step_ * xIndex;
	double y = yStart_ + step_ * yIndex;
	double z = zStart_ + step_ * zIndex;
	if (index.column() == 0) {
		return x;
	}
	else if (index.column() == 1) {
		return y;
	}
	else if (index.column() == 2) {
		return z;
	}
	else if (index.column() == 3) {
		return colorMap_->toColor(points_[zIndex][xIndex][yIndex]);
	}
	else if (index.column() == 4) {
		if ((xRate_ != Every && xIndex % xRate_) || (yRate_ != Every && yIndex % yRate_) || (zRate_ != Every && zIndex % zRate_) || checkXLimits(x) || checkYLimits(y) || checkZLimits(z))
		{
			return skippedPointSize_;
		}
		else
		{
			return pointSize_;
		}
	}
	else {
		return boost::any();
	}
	return boost::any();
}

bool VolumeData::checkXLimits(double x) const
{
	return x < limits_.xLeftLimit || x > limits_.xRightLimit;
}
bool VolumeData::checkYLimits(double y) const
{
	return y < limits_.yLeftLimit || y > limits_.yRightLimit;
}
bool VolumeData::checkZLimits(double z) const
{
	return z < limits_.zLeftLimit || z > limits_.zRightLimit;
}

VolumeData::~VolumeData()
{
	my_delete(colorMap_);
}
