#ifndef CHART_2D_PANEL_H
#define CHART_2D_PANEL_H

#include <Wt/Chart/WGridData>

#include "BaseChart.h"
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
	WGridData* gridData_;

	WText* pointValues_;
};

#endif

