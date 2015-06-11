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


SurfaceData2::SurfaceData2(FBspSurf& surf, unsigned nbXpts, unsigned nbYpts, WObject *parent)
: WStandardItemModel(nbXpts + 2, nbYpts + 2, parent),
xStart_(0), xEnd_(1), yStart_(0), yEnd_(1), xStep_((xEnd_ - xStart_) / (nbXpts)), yStep_((yEnd_ - yStart_) / (nbYpts))
{
	data_ = surf.ComputePoints(nbXpts, nbYpts);
}

boost::any SurfaceData2::data(const WModelIndex& index, int role) const
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
		return yStart_ + (index.column()-1)*yStep_;
	}
	if (index.column() == 0)
	{
		if (index.row() == 0)
		{
			return 0.0;
		}
		return xStart_ + (index.row()-1)*xStep_;
	}

	return data_.GetRow(index.row()-1)[index.column()-1];
}

VolumeData::VolumeData(Matrix3D<Point1D> points, double xStart, double yStart, double zStart, double step, Wt::WObject *parent)
: WStandardItemModel(points.size()*points.size()*points.size(), 4, parent), nbPts_(points.size()), points_(points), xStart_(xStart), yStart_(yStart), zStart_(zStart), step_(step), min_(INT_MAX), max_(INT_MIN)
{
	findMinAndMax();
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
	unsigned xIndex = ((index.row() % (nbPts_ * nbPts_)) / nbPts_);
	unsigned yIndex = (index.row() % nbPts_);
	unsigned zIndex = (index.row() / (nbPts_ * nbPts_));
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
	else {
		return boost::any();
	}
	return boost::any();
}

VolumeData::~VolumeData()
{
	delete colorMap_;
}
