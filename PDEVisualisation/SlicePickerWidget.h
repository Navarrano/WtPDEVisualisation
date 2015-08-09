#ifndef SLICE_PICKER_WIDGET_H
#define SLICE_PICKER_WIDGET_H

#include <Wt/WText>
#include <Wt/WContainerWidget>
#include <Wt/WComboBox>
#include <Wt/WDoubleSpinBox>
#include <Wt/WPushButton>

#include <boost/signal.hpp>

#include "structures.h"

using namespace Wt;

class SlicePickerWidget :
	public WContainerWidget
{
public:
	SlicePickerWidget(WContainerWidget* parent = 0);
	virtual ~SlicePickerWidget();
	void update(ChartData chartData, int precision = 2);
	boost::signal<void(Dimension, double)>& clicked();
private:
	WComboBox *comboBox_;
	WDoubleSpinBox *spinBox_;
	WPushButton *showButton_;
	WText *errorText_;
	double uMin_, vMin_, wMin_, uMax_, vMax_, wMax_;

	boost::signal<void(Dimension, double)> buttonClickedSignal;
	void updateSpinBox();
};

#endif
