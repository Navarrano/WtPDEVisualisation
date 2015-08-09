#include "Validators.h"

#include "PolynomialParser.h"

/////// DOUBLE PAIRED VALIDATOR

DoublePairedValidator::DoublePairedValidator(DoubleLineEdit* compareWith, Wt::WObject* parent) : Wt::WDoubleValidator(parent), compareWith_(compareWith) {}

Wt::WValidator::Result DoublePairedValidator::validate(const Wt::WString &input) const
{
	WValidator::Result r = WDoubleValidator::validate(input);
	if (r.state() != WValidator::Valid || !compareWith_)
		return r;
	if (compareWith_->validate() != WValidator::Valid)
		return WValidator::Result(WValidator::Invalid, "Minimum value is required");
	double val = compareWith_->value();
	if (std::stod(input.toUTF8()) <= val)
		return WValidator::Result(WValidator::Invalid, "Value has to be greater than minimum value");
	else
		return WValidator::Result(WValidator::Valid);
	
}

/////// POLYNOMIAL VALIDATOR

PolynomialValidator::PolynomialValidator(char* invalidVariable, Wt::WObject* parent) : WValidator(parent), invalidVariable_(invalidVariable) {}

Wt::WValidator::Result PolynomialValidator::validate(const Wt::WString &input) const
{
	WValidator::Result r = WValidator::validate(input);
	if (r.state() != WValidator::Valid)
		return r;
	PolynomialParser* p;
	if (invalidVariable_)
		p = new PolynomialParser(input.toUTF8(), *invalidVariable_);
	else
		p = new PolynomialParser(input.toUTF8());
	try
	{
		p->validate();
	}
	catch (std::invalid_argument& e)
	{
		return WValidator::Result(WValidator::Invalid, e.what());
	}
	return WValidator::Result(WValidator::Valid);
}
