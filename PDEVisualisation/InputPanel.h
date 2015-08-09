#ifndef INPUT_PANEL_H
#define INPUT_PANEL_H

#include "AbstractPanel.h"

#include <Wt/WPushButton>

#include <Wt/WLineEdit>
#include <Wt/WDoubleValidator>

#include "PolynomialParser.h"
#include "CustomGroupBoxes.h"

class InputPanel : public AbstractPanel
{
public:
	InputPanel(WContainerWidget* parent = 0);
	~InputPanel();
protected:
	void initComponents();
	void initEvents();
private:
	FBspVol preperVolData(SurfaceOrders orders, Limits limits);
	FBspVol calculatePDE(SurfaceOrders orders, Limits limits, std::array<PolynomialData, 6> boundaries, PolynomialData rhs);
	void fillVector(std::vector<int>& nums, int val, size_t start = 0, size_t lenght = 0);

	FBspVol pdeData_;

	OrdersGroupBox* ordersBox_;
	LimitsGroupBox* limitsBox_;
	BoundariesGroupBox* boundariesBox_;
	RhsGroupBox* rhsGroupBox_;
	SystemSettingsBox* settingsBox_;

	WPushButton* calculateButton_;
};



#endif


