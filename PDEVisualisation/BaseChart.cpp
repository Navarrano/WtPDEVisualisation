#include "BaseChart.h"

#include <Wt/WCssDecorationStyle>
#include <Wt/Chart/WStandardColorMap>
#include <Wt/Chart/WScatterData>

BaseChart::BaseChart(WContainerWidget *parent) : WCartesian3DChart(ScatterPlot, parent)
{
	this->setRenderOptions(Wt::WGLWidget::ClientSideRendering | Wt::WGLWidget::AntiAliasing);
	Wt::WCssDecorationStyle style;
	style.setBorder(Wt::WBorder(Wt::WBorder::Solid, Wt::WBorder::Medium,
		Wt::black));
	this->setDecorationStyle(style);

	this->setGridEnabled(Wt::Chart::XY_Plane, Wt::Chart::XAxis_3D, true);
	this->setGridEnabled(Wt::Chart::XY_Plane, Wt::Chart::YAxis_3D, true);
	this->setGridEnabled(Wt::Chart::XZ_Plane, Wt::Chart::XAxis_3D, true);
	this->setGridEnabled(Wt::Chart::XZ_Plane, Wt::Chart::ZAxis_3D, true);
	this->setGridEnabled(Wt::Chart::YZ_Plane, Wt::Chart::YAxis_3D, true);
	this->setGridEnabled(Wt::Chart::YZ_Plane, Wt::Chart::ZAxis_3D, true);

	this->axis(Wt::Chart::XAxis_3D).setTitle("u");
	this->axis(Wt::Chart::YAxis_3D).setTitle("v");
	this->axis(Wt::Chart::ZAxis_3D).setTitle("w");

	this->setIntersectionLinesEnabled(true);
	this->setIntersectionLinesColor(Wt::WColor(0, 255, 255));

	this->resize(607, 592);

	worldTransform_.lookAt(
		0.5, 0.5, 5, // camera position
		0.5, 0.5, 0.5,      // looking at
		0, 1, 0);        // 'up' vector

	this->setPerspectiveView();
}

void BaseChart::addDataset(WAbstractItemModel *model, bool meshEnabled)
{
	this->clearDatasets();

	WGridData* gridData_ = new WGridData(model);
	gridData_->setType(Wt::Chart::SurfaceSeries3D);
	gridData_->setSurfaceMeshEnabled(meshEnabled);
	Wt::Chart::WStandardColorMap *colorMap = new Wt::Chart::WStandardColorMap(gridData_->minimum(Wt::Chart::ZAxis_3D),
		gridData_->maximum(Wt::Chart::ZAxis_3D), true);
	gridData_->setColorMap(colorMap);
	this->addDataSeries(gridData_);

	gridData_->setColorMapVisible(true);
	gridData_->setColorMapSide(Left);
}

void BaseChart::changeAxisTitles(Dimension sliceDim)
{
	switch (sliceDim)
	{
	case U:
		this->axis(Wt::Chart::XAxis_3D).setTitle("v");
		this->axis(Wt::Chart::YAxis_3D).setTitle("w");
		this->axis(Wt::Chart::ZAxis_3D).setTitle("z");
		break;
	case V:
		this->axis(Wt::Chart::XAxis_3D).setTitle("u");
		this->axis(Wt::Chart::YAxis_3D).setTitle("w");
		this->axis(Wt::Chart::ZAxis_3D).setTitle("z");
		break;
	case W:
		this->axis(Wt::Chart::XAxis_3D).setTitle("u");
		this->axis(Wt::Chart::YAxis_3D).setTitle("v");
		this->axis(Wt::Chart::ZAxis_3D).setTitle("z");
		break;
	default:
		this->axis(Wt::Chart::XAxis_3D).setTitle("u");
		this->axis(Wt::Chart::YAxis_3D).setTitle("v");
		this->axis(Wt::Chart::ZAxis_3D).setTitle("w");
		break;
	}
}

void BaseChart::addScatterDataset(WAbstractItemModel *model, double min, double max)
{
	this->clearDatasets();

	WScatterData *scatter = new WScatterData(model);
	//scatter->setPointSize(7);

	WStandardColorMap *colorMap = new WStandardColorMap(min, max, true);
	scatter->setColorMap(colorMap);
	scatter->setColorColumn(3);
	scatter->setSizeColumn(4);

	this->addDataSeries(scatter);
	scatter->setColorMapVisible(true);
	scatter->setColorMapSide(Left);
}

void BaseChart::setTopView()
{
	WMatrix4x4 cameraMat = worldTransform_;
	cameraMat.translate(0.5, 0.5, 0.5);
	cameraMat.rotate(90.0, 1.0, 0.0, 0.0);
	cameraMat.scale(scale_);
	cameraMat.translate(-0.5, -0.5, -0.5);
	this->setCameraMatrix(cameraMat);
}

void BaseChart::setSideView()
{
	WMatrix4x4 cameraMat = worldTransform_;
	cameraMat.translate(0.5, 0.5, 0.5);
	cameraMat.rotate(90.0, 0.0, 1.0, 0.0);
	cameraMat.scale(scale_);
	cameraMat.translate(-0.5, -0.5, -0.5);
	this->setCameraMatrix(cameraMat);
}

void BaseChart::setFrontView()
{
	WMatrix4x4 cameraMat = worldTransform_;
	cameraMat.translate(0.5, 0.5, 0.5);
	cameraMat.rotate(0, 0.0, 1.0, 0.0);
	cameraMat.scale(scale_);
	cameraMat.translate(-0.5, -0.5, -0.5);
	this->setCameraMatrix(cameraMat);
}

void BaseChart::setPerspectiveView()
{
	WMatrix4x4 cameraMat = worldTransform_;
	cameraMat.translate(0.5, 0.5, 0.5);
	cameraMat.rotate(45.0, 0.0, 1.0, 0.0);
	cameraMat.rotate(20.0, 1.0, 0.0, 1.0);
	cameraMat.scale(scale_);
	cameraMat.translate(-0.5, -0.5, -0.5);
	this->setCameraMatrix(cameraMat);
}

void BaseChart::changeView(ViewType viewType)
{
	switch (viewType)
	{
	case ViewType::TOP:
		setTopView();
		break;
	case ViewType::SIDE:
		setSideView();
		break;
	case ViewType::FRONT:
		setFrontView();
		break;
	case ViewType::PERSPECTIVE:
		setPerspectiveView();
		break;
	default:
		setPerspectiveView();
		break;
	}
}

void BaseChart::clearDatasets()
{
	const std::vector<WAbstractDataSeries3D*> &vec = this->dataSeries();
	std::for_each(vec.begin(), vec.end(), [&](WAbstractDataSeries3D* series) {removeDataSeries(series); });
}

BaseChart::~BaseChart()
{
}
