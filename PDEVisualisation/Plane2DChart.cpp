#include "Plane2DChart.h"
#include "Models.h"

#include <Wt/WCssDecorationStyle>
#include <Wt/Chart/WGridData>
#include <Wt/Chart/WStandardColorMap>


Plane2DChart::Plane2DChart(WContainerWidget *parent) : BaseChart(parent)
{
	this->addStyleClass("disabled");
	changeView(ViewType::TOP);
}

Plane2DChart::~Plane2DChart()
{
}
