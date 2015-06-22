#include "InputPanel.h"

#include <Wt/WText>
#include <Wt/WGroupBox>
#include <Wt/WLabel>

InputPanel::InputPanel(WContainerWidget* parent) : AbstractPanel("Input data", parent)
{
	initComponents();

	std::ifstream ifs("F:/LaplaceSoln.dat");
	ifs >> pdeData_;
}

void InputPanel::initComponents()
{
	WGroupBox* pdeBox = new WGroupBox("Partial derivative equation", root());
	WGroupBox *groupBox = new WGroupBox("Charts settings", root());
	groupBox->setMargin(20, Bottom);
	WLabel *densityLabel = new WLabel("Density: ", groupBox);
	WComboBox *densityBox = new WComboBox(groupBox);
	for (int i = 0; i < stepsArraySize; i++)
	{
		densityBox->addItem(std::to_string(possibleDensities[i]));
	}
	densityBox->setCurrentIndex(3);

	new WText("<hr/>", root());
	calculateButton_ = new WPushButton("Calculate", root());
	calculateButton_->clicked().connect(std::bind([=](){
		ChartData chartData;
		unsigned density = possibleDensities[densityBox->currentIndex()];
		chartData.xStart = pdeData_.GetLeftLimitU();
		chartData.xEnd = pdeData_.GetRightLimitU();
		chartData.yStart = pdeData_.GetLeftLimitV();
		chartData.yEnd = pdeData_.GetRightLimitV();
		chartData.zStart = pdeData_.GetLeftLimitW();
		chartData.zEnd = pdeData_.GetRightLimitW();
		chartData.xSize = (chartData.xEnd - chartData.xStart) * density + 1;
		chartData.ySize = (chartData.yEnd - chartData.yStart) * density + 1;
		chartData.zSize = (chartData.zEnd - chartData.zStart) * density + 1;
		chartData.totalNbPts = chartData.xSize * chartData.ySize * chartData.zSize;
		chartData.points = pdeData_.ComputePoints(chartData.xSize - 1, chartData.ySize - 1, chartData.zSize - 1);
		chartData.step = possibleStepSizes[densityBox->currentIndex()];
		activate(chartData);
		//activate(pdeData_, chartData);
	}));
}

InputPanel::~InputPanel()
{
	delete calculateButton_;
}
