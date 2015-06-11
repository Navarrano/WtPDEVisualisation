#include "InputPanel.h"

#include <Wt/WText>

InputPanel::InputPanel(WContainerWidget* parent) : AbstractPanel("Input data", parent)
{
	initComponents();

	std::ifstream ifs("F:/LaplaceSoln.dat");
	ifs >> pdeData_;
}

void InputPanel::initComponents()
{
	root()->addWidget(new WText("<h2>Input panel</h2>", root()));

	calculateButton_ = new WPushButton("Calculate", root());
	unsigned size = 100;
	calculateButton_->clicked().connect(std::bind([=](){
		ChartData chartData;
		chartData.points = pdeData_.ComputePoints(size, size, size);
		chartData.nbPts = size;
		chartData.xStart = 0.0;
		chartData.xEnd = 1.0;
		chartData.yStart = 0.0;
		chartData.yEnd = 1.0;
		chartData.zStart = 0.0;
		chartData.zEnd = 1.0;
		chartData.nbPts = size;
		chartData.step = 0.01;
		activate(chartData);
	}));
}

InputPanel::~InputPanel()
{
	delete calculateButton_;
}
