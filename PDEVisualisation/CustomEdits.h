#ifndef CUSTOM_EDITS_H
#define CUSTOM_EDITS_H

#include <Wt/WLineEdit>
#include <boost/signal.hpp>

#include "structures.h"
#include "PolynomialParser.h"

/////// DOUBLE LINE EDIT

class DoublePairedValidator;

class DoubleLineEdit : public Wt::WLineEdit
{
public:
	DoubleLineEdit(DoubleLineEdit* compareWith = 0, Wt::WContainerWidget* parent = 0);
	virtual ~DoubleLineEdit() {};
	const double value();
	void setValue(const double val);
	boost::signal<void(DoubleLineEdit*)>& changed() { return signal_; };
private:
	boost::signal<void(DoubleLineEdit*)> signal_;
};

/////// POLYNOMIAL LINE EDIT

class PolynomialLineEdit : public Wt::WLineEdit
{
public:
	PolynomialLineEdit(char* invalidVar = 0, Wt::WContainerWidget* parent = 0);
	virtual ~PolynomialLineEdit() {};
	PolynomialData value();
private:
	char* invalidVariable_;
};

#endif