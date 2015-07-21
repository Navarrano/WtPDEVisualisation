#include "Chart2DPanel.h"
#include "Models.h"
#include "ViewButtonsWidget.h"

#include <Wt/WText>
#include <Wt/WVBoxLayout>
#include <Wt/WHBoxLayout>
#include <Wt/WLineEdit>
#include <string>

std::string doubleToString(const double val, short precision = 2)
{
	std::stringstream str;
	str << std::fixed << std::setprecision(precision) << val;
	return str.str();
}

Chart2DPanel::Chart2DPanel(WContainerWidget *parent) : AbstractPanel("Slice Chart", parent), gridData_(0)
{
	this->collapse();
	initComponents();
}

void Chart2DPanel::initComponents()
{
	WVBoxLayout *vbox = new WVBoxLayout();
	WHBoxLayout *hbox = new WHBoxLayout();
	WContainerWidget *firstChart = new WContainerWidget();
	WContainerWidget *secondChart = new WContainerWidget();
	WContainerWidget *comboBoxContainer = new WContainerWidget();

	planeChart_ = new Plane2DChart(firstChart);
	pointValues_ = new WText("", firstChart);
	surfaceChart_ = new BaseChart(secondChart);

	ViewButtonsWidget *viewButtons = new ViewButtonsWidget(secondChart); 
	viewButtons->clicked().connect(boost::bind(&BaseChart::changeView, surfaceChart_, _1));

	hbox->addWidget(firstChart);
	hbox->addWidget(secondChart);
	slicePicker_ = new SlicePickerWidget(comboBoxContainer);
	slicePicker_->clicked().connect(boost::bind(&Chart2DPanel::redrawCharts, this, _1, _2));

	comboBoxContainer->addWidget(slicePicker_);

	vbox->addLayout(hbox, 1);
	vbox->addWidget(comboBoxContainer);
	root()->setLayout(vbox);

	planeChart_->clicked().connect(std::bind([&](WMouseEvent& e)
	{
		int x = e.widget().x;
		int y = e.widget().y;
		if (gridData_)
		{
			std::vector<WSurfaceSelection> surfaces = gridData_->pickSurface(x, y);
			if (surfaces.size())
			{
				WSurfaceSelection surf = surfaces[0];
				WString str("f({1},{2}) = {3}");
				pointValues_->setText(str.arg(doubleToString(surf.x)).arg(doubleToString(surf.y)).arg(doubleToString(surf.z)));
			}
		}
	}, std::placeholders::_1));
}

void Chart2DPanel::update(ChartData& chartData)
{
	chartData_ = chartData;
	slicePicker_->update(chartData);
	redrawCharts(Dimension::U, chartData.xStart);
	AbstractPanel::update(chartData);
}

unsigned Chart2DPanel::valueToIndex(double min, double step, double value)
{
	return (value - min) / step;
}

void Chart2DPanel::redrawCharts(Dimension dim, double value)
{
	Matrix<Point1D> points;
	unsigned index;
	double xStart, yStart;
	switch (dim)
	{
	case U:
		index = valueToIndex(chartData_.xStart, chartData_.step, value);
		xStart = chartData_.yStart;
		yStart = chartData_.zStart;
		points = chartData_.points.GetVW(index);
		break;
	case V:
		index = valueToIndex(chartData_.yStart, chartData_.step, value);
		xStart = chartData_.xStart;
		yStart = chartData_.zStart;
		points = chartData_.points.GetUW(index);
		break;
	case W:
		index = valueToIndex(chartData_.zStart, chartData_.step, value);
		xStart = chartData_.xStart;
		yStart = chartData_.yStart;
		points = chartData_.points.GetUV(index);
		break;
	default:
		break;
	}
	planeChart_->changeAxisTitles(dim);
	gridData_ = planeChart_->addDataset(new SurfaceData(points, xStart, yStart, chartData_.step, planeChart_), false);
	
	surfaceChart_->changeAxisTitles(dim);
	surfaceChart_->addDataset(new SurfaceData(points, xStart, yStart, chartData_.step, surfaceChart_));
}

Chart2DPanel::~Chart2DPanel()
{
	delete planeChart_;
	delete surfaceChart_;
	delete slicePicker_;
	delete pointValues_;
}
