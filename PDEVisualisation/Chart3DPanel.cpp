#include "Chart3DPanel.h"

#include <Wt/WText>
#include <Wt/WHBoxLayout>
#include <Wt/WLabel>

#include "Models.h"
#include "ViewButtonsWidget.h"



Chart3DPanel::Chart3DPanel(WContainerWidget *parent) : AbstractPanel("3D Volume Chart", parent)
{
	this->collapse();

	initComponents();
}

void Chart3DPanel::initComponents()
{
	WHBoxLayout *hbox = new WHBoxLayout();
	WContainerWidget *chartContainer = new WContainerWidget();
	WContainerWidget *settingsContainer = new WContainerWidget();

	volumeChart_ = new BaseChart(chartContainer);
	volumeChart_->resize(846, 846);
	//volumeChart_->setBackground(WColor(235, 245, 223));
	//volumeChart_->setBackground(WColor(135, 135, 135));
	settingsBox_ = new SettingsBox(settingsContainer);

	ViewButtonsWidget *viewButtons = new ViewButtonsWidget(chartContainer);
	viewButtons->clicked().connect(boost::bind(&BaseChart::changeView, volumeChart_, _1));

	hbox->addWidget(chartContainer);
	hbox->addWidget(settingsContainer, 1);

	root()->setLayout(hbox);
}

void Chart3DPanel::update(ChartData& chartData)
{
	currentChartData_ = chartData;
	VolumeData* volumeData = new VolumeData(chartData.points, chartData.xStart, chartData.yStart, chartData.zStart, chartData.step, volumeChart_);
	
	volumeChart_->clearDatasets();
	volumeChart_->addScatterDataset(volumeData, volumeData->getMin(), volumeData->getMax());

	settingsBox_->update(chartData);
	AbstractPanel::update(chartData);
}

Chart3DPanel::~Chart3DPanel()
{
	delete volumeChart_;
	delete settingsBox_;
}

SettingsBox::SettingsBox(WContainerWidget *parent) : WGroupBox("Volume chart settings", parent)
{
	initComponenets();
	initEvents();
}

void SettingsBox::initComponenets()
{
	WLabel *densityLabel = new WLabel("Density: ", this);
	densityBox_ = new WComboBox(this);
	densityBox_->setInline(true);
	densityBox_->setDisabled(true);
	densityLabel->setBuddy(densityBox_);
}

void SettingsBox::initEvents()
{
	densityBox_->changed().connect(std::bind([=](){
		densityChangedSignal_(densityBox_->currentIndex());
	}));
}

void SettingsBox::update(ChartData& chartData)
{
	activateControls();
	updateDensityBox(chartData.step);
}

void SettingsBox::activateControls()
{
	densityBox_->setDisabled(false);
}

void SettingsBox::updateDensityBox(double maxDensity)
{
	for (int i = 0; i < stepsArraySize; i++)
	{
		if (possibleStepSizes[i] >= maxDensity)
		{
			densityBox_->addItem(std::to_string(possibleDensities[i]));		
		}	
	}
	densityBox_->setCurrentIndex(densityBox_->count());
}

SettingsBox::~SettingsBox() 
{
	delete densityBox_;
}
