#ifndef VIEW_BUTTONS_WIDGET_H
#define VIEW_BUTTONS_WIDGET_H

#include "BaseChart.h"

#include <Wt/WContainerWidget>
#include <Wt/WPushButton>
#include <boost/signal.hpp>

typedef boost::signal<void(ViewType)> ViewSignal;

class ViewButtonsWidget :
	public WContainerWidget
{
public:
	ViewButtonsWidget(WContainerWidget* parent = 0);
	ViewSignal& clicked() { return viewChangeSignal_; };
	virtual ~ViewButtonsWidget();
private:
	WPushButton* perspectiveView_;
	WPushButton* topView_;
	WPushButton* sideView_;
	WPushButton* frontView_;
	ViewSignal viewChangeSignal_;
};

#endif

