#include "CustomEdits.h"

#include "Validators.h"
#include "global_methods.h"

/////// DOUBLE LINE EDIT

DoubleLineEdit::DoubleLineEdit(DoubleLineEdit* compareWith, Wt::WContainerWidget* parent)
{
	Wt::WValidator* validator = new DoublePairedValidator(compareWith, this);
	validator->setMandatory(true);
	this->setValidator(validator);
	WLineEdit::blurred().connect(std::bind([=]() { signal_(this); }));
}

const double DoubleLineEdit::value()
{
	if (validate() == Wt::WValidator::Valid)
	{
		std::string str = (this->valueText()).toUTF8();
		return std::stod(str);
	}
	throw std::invalid_argument("Can't convert field value to double");
}

void DoubleLineEdit::setValue(const double val)
{
	this->setValueText(toString(val));
}

/////// POLYNOMIAL LINE EDIT

PolynomialLineEdit::PolynomialLineEdit(char* invalidVar, Wt::WContainerWidget* parent) : Wt::WLineEdit(parent), invalidVariable_(invalidVar)
{
	Wt::WValidator* validator = new PolynomialValidator(invalidVariable_, this);
	validator->setMandatory(true);
	this->setValidator(validator);
	//WLineEdit::blurred().connect(std::bind([=]() { this->validate(); }));
}

PolynomialData PolynomialLineEdit::value()
{
	PolynomialParser* p;
	if (invalidVariable_)
		p = new PolynomialParser((this->valueText()).toUTF8(), *invalidVariable_);
	else
		p = new PolynomialParser((this->valueText()).toUTF8());
	return p->parse();
}