#ifndef CUSTOM_GROUP_BOXES_H
#define CUSTOM_GROUP_BOXES_H

#include "structures.h"
#include "global_methods.h"
#include "CustomEdits.h"

#include <Wt/WGroupBox>
#include <Wt/WComboBox>
#include <Wt/WText>
#include <boost/signal.hpp>


//////// ORDERS GROUP BOX ////////

class OrdersGroupBox : public Wt::WGroupBox
{
public:
	OrdersGroupBox(Wt::WContainerWidget* parent = 0);
	virtual ~OrdersGroupBox();
	SurfaceOrders getOrders();
private:
	void initComponents();
	Wt::WComboBox* initComboBox();

	Wt::WComboBox* xOrder_;
	Wt::WComboBox* yOrder_;
	Wt::WComboBox* zOrder_;
};

//////// LIMITS GROUP BOX ////////

class LimitsGroupBox : public Wt::WGroupBox
{
public:
	LimitsGroupBox(Wt::WContainerWidget* parent = 0);
	virtual ~LimitsGroupBox();
	boost::signal<void()>& valid() { return validSignal_; };
	boost::signal<void()>& invalid() { return invalidSignal_; };
	Limits getLimits();
private:
	void initComponents();
	void initEvents();
	bool areAllValid();

	DoubleLineEdit* inputs_[6];

	std::array<bool, 6> valid_;
	boost::signal<void()> validSignal_;
	boost::signal<void()> invalidSignal_;
};


//////// BOUNDARIES GROUP BOX ////////

class BoundariesGroupBox : public Wt::WGroupBox
{
public:
	BoundariesGroupBox(Wt::WContainerWidget* parent = 0);
	virtual ~BoundariesGroupBox();
	void enableAll();
	void disableAll();
	void update(Limits limits);
	bool enabled() { return enabled_; };
	std::array<PolynomialData, 6> boundaries();
	boost::signal<void()>& valid() { return validSignal_; };
	boost::signal<void()>& invalid() { return invalidSignal_; };
private:
	void initComponents();
	void initEvents();
	void initEditEvent(Wt::WLineEdit* edit, size_t idx);
	bool areAllValid();
	Wt::WString prepareLabel(std::string u, std::string v, std::string w);

	PolynomialLineEdit* xLeft_;
	Wt::WText* xLeftLabel_;
	PolynomialLineEdit* xRight_;
	Wt::WText* xRightLabel_;

	PolynomialLineEdit* yLeft_;
	Wt::WText* yLeftLabel_;
	PolynomialLineEdit* yRight_;
	Wt::WText* yRightLabel_;

	PolynomialLineEdit* zLeft_;
	Wt::WText* zLeftLabel_;
	PolynomialLineEdit* zRight_;
	Wt::WText* zRightLabel_;

	bool enabled_;
	std::array<bool, 6> valid_;
	boost::signal<void()> validSignal_;
	boost::signal<void()> invalidSignal_;
};

//////// RHS GROUP BOX ////////

class RhsGroupBox : public Wt::WGroupBox
{
public:
	RhsGroupBox(Wt::WContainerWidget* parent = 0);
	virtual ~RhsGroupBox();
	PolynomialData rhs() { return rhs_->value(); };
	void enable();
	void disable() { rhs_->setDisabled(true); };
	boost::signal<void()>& valid() { return validSignal_; };
	boost::signal<void()>& invalid() { return invalidSignal_; };
private:
	void initComponents();
	void initEvents();

	PolynomialLineEdit* rhs_;
	boost::signal<void()> validSignal_;
	boost::signal<void()> invalidSignal_;
};

//////// SYSTEM SETTINGS BOX ////////
class SystemSettingsBox : public Wt::WGroupBox
{
public:
	SystemSettingsBox(Wt::WContainerWidget* parent = 0);
	virtual ~SystemSettingsBox();
	unsigned density();
	double stepSize();
private:
	void initComponents();

	Wt::WComboBox* density_;
};

#endif