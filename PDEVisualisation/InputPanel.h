#pragma once

#include "AbstractPanel.h"

#include <Wt/WPushButton>
#include <Wt/WGroupBox>
#include <Wt/WLineEdit>
#include <Wt/WDoubleValidator>

#include "PolynomialParser.h"

class LimitsGroupBox;
class BoundariesGroupBox;
class DoubleLineEdit;
class ResultBox;
class PolynomialLineEdit;

class InputPanel : public AbstractPanel
{
public:
	InputPanel(WContainerWidget* parent = 0);
	~InputPanel();
protected:
	void initComponents();
	void initEvents();
private:
	FBspVol& preperVolData(SurfaceOrders orders, Limits limits);
	FBspVol& calculatePDE(SurfaceOrders orders, Limits limits, std::array<PolynomialData, 6> boundaries, PolynomialData rhs);
	WComboBox* initComboBox();

	WPushButton* calculateButton_;
	FBspVol pdeData_;
	LimitsGroupBox* limitsBox_;
	BoundariesGroupBox* boundariesBox_;

	WComboBox* xOrder_;
	WComboBox* yOrder_;
	WComboBox* zOrder_;
	PolynomialLineEdit* rhs_;
	WComboBox* density_;
};

///// CUSTOM VALIDATORS

class CustomValidator : public WDoubleValidator
{
public:
	CustomValidator(DoubleLineEdit* compareWith, WObject* parent = 0);
	virtual ~CustomValidator() {};
	virtual WValidator::Result validate(const WString &input) const;
private:
	DoubleLineEdit* compareWith_;
};

class PolynomialValidator : public WValidator
{
public:
	PolynomialValidator(WObject* parent = 0);
	PolynomialValidator(const char invalidVariable, WObject* parent = 0);
	virtual ~PolynomialValidator() {};
	virtual WValidator::Result validate(const WString &input) const;
private:
	bool twoVar_;
	char invalidVariable_;
};

//// DOUBLE LINE EDIT

class DoubleLineEdit : public WLineEdit
{
public:
	DoubleLineEdit(WContainerWidget* parent = 0);
	DoubleLineEdit(DoubleLineEdit* compareWith, WContainerWidget* parent = 0);
	virtual ~DoubleLineEdit();
	const double value();
	void setValue(const double val);
	boost::signal<void(DoubleLineEdit*)>& changed() { return signal_; };
private:
	boost::signal<void(DoubleLineEdit*)> signal_;
};

/////// POLYNOMIAL LINE EDIT

class PolynomialLineEdit : public WLineEdit
{
public:
	PolynomialLineEdit(WContainerWidget* parent = 0);
	PolynomialLineEdit(const char invalidVariable, WContainerWidget* parent = 0);
	virtual ~PolynomialLineEdit() {};
	PolynomialData value();
	boost::signal<void(PolynomialLineEdit*)>& changed() { return signal_; };
private:
	boost::signal<void(PolynomialLineEdit*)> signal_;
	bool twoVar_;
	char invalidVariable_;
};

//// LIMITS GROUP BOX

class LimitsGroupBox : public WGroupBox
{
public:
	LimitsGroupBox(WContainerWidget* parent = 0);
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


////// BOUNDARIES GROUP BOX

class BoundariesGroupBox : public WGroupBox
{
public:
	BoundariesGroupBox(WContainerWidget* parent = 0);
	virtual ~BoundariesGroupBox();
	void enableAll();
	void disableAll();
	void update(Limits limits);
	bool enabled() { return enabled_; };
	std::array<PolynomialData, 6> boundaries();
private:
	void initComponents();
	void initEvents();
	WString prepareLabel(std::string u, std::string v, std::string w);

	PolynomialLineEdit* xLeft_;
	WText* xLeftLabel_;
	PolynomialLineEdit* xRight_;
	WText* xRightLabel_;

	PolynomialLineEdit* yLeft_;
	WText* yLeftLabel_;
	PolynomialLineEdit* yRight_;
	WText* yRightLabel_;

	PolynomialLineEdit* zLeft_;
	WText* zLeftLabel_;
	PolynomialLineEdit* zRight_;
	WText* zRightLabel_;

	bool enabled_;
};


/// RESULT BOX

class ResultBox : public WGroupBox
{
public:
	ResultBox(WContainerWidget* parent = 0);
	virtual ~ResultBox();
private:
	WLineEdit* lhs_;
};




