#pragma once

#include "Plane2DChart.h"
#include "AbstractPanel.h"
#include "SlicePickerWidget.h"

class Chart2DPanel :
	public AbstractPanel
{
public:
	Chart2DPanel(WContainerWidget *parent = 0);
	virtual ~Chart2DPanel();
	virtual void update(ChartData& chartData);
private:
	void initComponents();
	void redrawCharts(Dimension dim, double value);
	unsigned valueToIndex(double min, double step, double value);

	BaseChart* planeChart_;
	BaseChart* surfaceChart_;
	SlicePickerWidget* slicePicker_;
	ChartData chartData_;
};

