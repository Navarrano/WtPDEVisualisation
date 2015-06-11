#pragma once

#include <Wt/WContainerWidget>
#include <Wt/Chart/WCartesian3DChart>
#include <Wt/WAbstractItemModel>
#include <Wt/Chart/WGridData>

#include "headers.h"

using namespace Wt;
using namespace Wt::Chart;

enum ViewType { TOP, SIDE, FRONT, PERSPECTIVE };

class BaseChart :
	public WCartesian3DChart
{
public:
	BaseChart(WContainerWidget* parent = 0);
	virtual ~BaseChart();
	void addDataset(WAbstractItemModel* model, bool meshEnabled = true);
	void addScatterDataset(WAbstractItemModel *model, double min, double max);
	void changeView(ViewType viewType = ViewType::PERSPECTIVE);
	void changeAxisTitles(Dimension sliceDim);
	void clearDatasets();

private:
	void setTopView();
	void setSideView();
	void setPerspectiveView();
	void setFrontView();

	const double scale_ = 2;
	WMatrix4x4 worldTransform_;
};

