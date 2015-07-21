#pragma once

#include <Wt/WContainerWidget>
#include <Wt/Chart/WCartesian3DChart>
#include <Wt/WAbstractItemModel>
#include <Wt/Chart/WGridData>
#include <Wt/Chart/WScatterData>

#include "headers.h"
#include "Models.h"

using namespace Wt;
using namespace Wt::Chart;

enum ViewType { TOP, SIDE, FRONT, PERSPECTIVE };

class BaseChart :
	public WCartesian3DChart
{
public:
	BaseChart(WContainerWidget* parent = 0);
	virtual ~BaseChart();
	virtual WGridData* addDataset(WAbstractItemModel* model, bool meshEnabled = true);
	WScatterData* addVolumeDataset(VolumeData *model);
	void changeView(ViewType viewType = ViewType::PERSPECTIVE);
	void changeAxisTitles(Dimension sliceDim);
	void clearDatasets();
	void toggleColorMap(bool show);
private:
	void setTopView();
	void setSideView();
	void setPerspectiveView();
	void setFrontView();

	const double scale_ = 2;
	WMatrix4x4 worldTransform_;
};

class Plane2DChart :
	public BaseChart
{
public:
	Plane2DChart(WContainerWidget *parent = 0);
	virtual ~Plane2DChart();
};

