#pragma once

#include <Wt/WGroupBox>
#include <Wt/WComboBox>

#include <boost/signals.hpp>

#include "BaseChart.h"
#include "AbstractPanel.h"

class SettingsBox;

class Chart3DPanel :
	public AbstractPanel
{
public:
	Chart3DPanel(WContainerWidget *parent = 0);
	virtual void update(ChartData& chartData);
	virtual ~Chart3DPanel();
private:
	void initComponents();

	BaseChart* volumeChart_;
	SettingsBox* settingsBox_;
	ChartData currentChartData_;
};

class SettingsBox : public WGroupBox
{
public:
	SettingsBox(WContainerWidget *parent = 0);
	void update(ChartData& chartData);
	boost::signal<void(int)>& densityChanged() { return densityChangedSignal_; };
	virtual ~SettingsBox();
private:
	void initComponenets();
	void initEvents();
	void activateControls();
	void updateDensityBox(double maxDensity);

	WComboBox *densityBox_;
	boost::signal<void(int)> densityChangedSignal_;
};

