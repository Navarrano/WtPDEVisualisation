#pragma once

#include "AbstractPanel.h"

#include <Wt/WPushButton>

class InputPanel : public AbstractPanel
{
public:
	InputPanel(WContainerWidget* parent = 0);
	~InputPanel();
protected:
	void initComponents();
private:
	WPushButton* calculateButton_;
	FBspVol pdeData_;
};

