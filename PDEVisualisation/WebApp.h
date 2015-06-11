#pragma once

#include "InputPanel.h"
#include "Chart2DPanel.h"
#include "Chart3DPanel.h"


#include <Wt/WApplication>
#include <Wt/WContainerWidget>
#include <Wt/WEnvironment>

using namespace Wt;

class WebApp : public WApplication
{
public:
	WebApp(const WEnvironment& env);
	~WebApp();
private:
	WContainerWidget* wrapper_;
	AbstractPanel* inputPanel_;
	AbstractPanel* slicePanel_;
	AbstractPanel* volumePanel_;

	void initComponents();
};

