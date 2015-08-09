#ifndef VALIDATORS_H
#define VALIDATORS_H

#include <Wt/WDoubleValidator>
#include "CustomEdits.h"

class DoubleLineEdit;

class DoublePairedValidator : public Wt::WDoubleValidator
{
public:
	DoublePairedValidator(DoubleLineEdit* compareWith = 0, Wt::WObject* parent = 0);
	virtual ~DoublePairedValidator() {};
	virtual WValidator::Result validate(const Wt::WString &input) const;
private:
	DoubleLineEdit* compareWith_;
};

class PolynomialValidator : public Wt::WValidator
{
public:
	PolynomialValidator(char* invalidVariable = 0, Wt::WObject* parent = 0);
	virtual ~PolynomialValidator() { if (invalidVariable_) delete invalidVariable_;  };
	virtual Wt::WValidator::Result validate(const Wt::WString &input) const;
private:
	char* invalidVariable_;
};

#endif
