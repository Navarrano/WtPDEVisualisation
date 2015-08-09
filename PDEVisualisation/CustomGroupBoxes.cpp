#include "CustomGroupBoxes.h"

#include <Wt/WTemplate>

/////// ORDERS GROUP BOX

OrdersGroupBox::OrdersGroupBox(WContainerWidget* parent) : WGroupBox("Orders of surfaces", parent)
{
	setStyleClass("child-legend");
	initComponents();
}

void OrdersGroupBox::initComponents()
{

	Wt::WTemplate* t = new Wt::WTemplate(Wt::WString::tr("surfaces-form"), this);
	t->addFunction("id", &Wt::WTemplate::Functions::id);
	xOrder_ = initComboBox();
	t->bindWidget("xSurface", xOrder_);
	yOrder_ = initComboBox();
	t->bindWidget("ySurface", yOrder_);
	zOrder_ = initComboBox();
	t->bindWidget("zSurface", zOrder_);
}

Wt::WComboBox* OrdersGroupBox::initComboBox()
{
	Wt::WComboBox* cb = new Wt::WComboBox();
	for (int i = 2; i <= 6; i++)
	{
		cb->addItem(std::to_string(i));
	}
	cb->setCurrentIndex(1);
	return cb;
}

SurfaceOrders OrdersGroupBox::getOrders()
{
	return SurfaceOrders(xOrder_->currentIndex() + 2, yOrder_->currentIndex() + 2, zOrder_->currentIndex() + 2);
}

OrdersGroupBox::~OrdersGroupBox()
{
	my_delete(xOrder_);
	my_delete(yOrder_);
	my_delete(zOrder_);
}

/////// LIMITS GROUP BOX

LimitsGroupBox::LimitsGroupBox(WContainerWidget* parent) : WGroupBox("Limits", parent)
{
	this->addStyleClass("child-legend limits");

	valid_.fill(true);

	initComponents();
	initEvents();
}

void LimitsGroupBox::initComponents()
{
	Wt::WTemplate* t = new Wt::WTemplate(Wt::WString::tr("limits-form"), this);
	t->addFunction("id", &Wt::WTemplate::Functions::id);

	inputs_[0] = new DoubleLineEdit();
	inputs_[0]->setValue(0.0);
	inputs_[1] = new DoubleLineEdit();
	inputs_[1]->setValue(0.0);
	inputs_[2] = new DoubleLineEdit();
	inputs_[2]->setValue(0.0);
	inputs_[3] = new DoubleLineEdit(inputs_[0]);
	inputs_[3]->setValue(1.0);
	inputs_[4] = new DoubleLineEdit(inputs_[1]);
	inputs_[4]->setValue(1.0);
	inputs_[5] = new DoubleLineEdit(inputs_[2]);
	inputs_[5]->setValue(1.0);
	t->bindWidget("leftX", inputs_[0]);
	t->bindWidget("leftY", inputs_[1]);
	t->bindWidget("leftZ", inputs_[2]);
	t->bindWidget("rightX", inputs_[3]);
	t->bindWidget("rightY", inputs_[4]);
	t->bindWidget("rightZ", inputs_[5]);
}

Limits LimitsGroupBox::getLimits()
{
	return Limits(inputs_[0]->value(), inputs_[3]->value(), inputs_[1]->value(), inputs_[4]->value(), inputs_[2]->value(), inputs_[5]->value());
}

bool LimitsGroupBox::areAllValid()
{
	for (int i = 0; i < 6; ++i)
	{
		if (!valid_[i])
			return false;
	}
	return true;
}

void LimitsGroupBox::initEvents()
{
	for (int i = 0; i < 6; i++)
	{
		inputs_[i]->changed().connect(std::bind([=](){
			if (inputs_[i]->validate() == Wt::WValidator::Valid)
			{
				valid_[i] = true;
				if (i <= 2)
				{
					if (inputs_[i + 3]->validate() == Wt::WValidator::Valid)
						valid_[i + 3] = true;
					else
						valid_[i + 3] = false;
				}
			}
			else
				valid_[i] = false;
			if (areAllValid())
				validSignal_();
			else
				invalidSignal_();
		}));
	}
}

LimitsGroupBox::~LimitsGroupBox()
{
	for (int i = 0; i < 6; i++)
		my_delete(inputs_[i]);
}

////// BOUNDARIES GROUP BOX

BoundariesGroupBox::BoundariesGroupBox(WContainerWidget* parent) : WGroupBox("Boundaries", parent), enabled_(true)
{
	this->addStyleClass("child-legend");

	valid_.fill(false);

	initComponents();
	initEvents();
}

void BoundariesGroupBox::initComponents()
{
	Wt::WTemplate* t = new Wt::WTemplate(Wt::WString::tr("boundaries-form"), this);
	t->addFunction("id", &Wt::WTemplate::Functions::id);

	xLeft_ = new PolynomialLineEdit(new char('u'));
	xLeftLabel_ = new Wt::WText("F(0.00 ,v, w) =");
	xRight_ = new PolynomialLineEdit(new char('u'));
	xRightLabel_ = new Wt::WText("F(1.00 ,v ,w) =");

	yLeft_ = new PolynomialLineEdit(new char('v'));
	yLeftLabel_ = new Wt::WText("F(u, 0.00, w) =");
	yRight_ = new PolynomialLineEdit(new char('v'));
	yRightLabel_ = new Wt::WText("F(u, 1.00, w) =");

	zLeft_ = new PolynomialLineEdit(new char('w'));
	zLeftLabel_ = new Wt::WText("F(u, v, 0.00) =");
	zRight_ = new PolynomialLineEdit(new char('w'));
	zRightLabel_ = new Wt::WText("F(u, v, 1.00) =");

	t->bindWidget("xLeft", xLeft_);
	t->bindWidget("xLeftLabel", xLeftLabel_);
	t->bindWidget("xRight", xRight_);
	t->bindWidget("xRightLabel", xRightLabel_);

	t->bindWidget("yLeft", yLeft_);
	t->bindWidget("yLeftLabel", yLeftLabel_);
	t->bindWidget("yRight", yRight_);
	t->bindWidget("yRightLabel", yRightLabel_);

	t->bindWidget("zLeft", zLeft_);
	t->bindWidget("zLeftLabel", zLeftLabel_);
	t->bindWidget("zRight", zRight_);
	t->bindWidget("zRightLabel", zRightLabel_);
}

void BoundariesGroupBox::initEvents()
{
	initEditEvent(xLeft_, 0);
	initEditEvent(xRight_, 1);
	initEditEvent(yLeft_, 2);
	initEditEvent(yRight_, 3);
	initEditEvent(zLeft_, 4);
	initEditEvent(zRight_, 5);
}

void BoundariesGroupBox::initEditEvent(Wt::WLineEdit* edit, size_t idx)
{
	edit->blurred().connect(std::bind([=]() {
		if (edit->validate() == Wt::WValidator::Valid)
		{
			valid_[idx] = true;
		}
		else
		{
			valid_[idx] = false;
		}
		if (areAllValid())
			validSignal_();
		else
			invalidSignal_();
	}));
}

bool BoundariesGroupBox::areAllValid()
{
	for (int i = 0; i < 6; ++i)
	{
		if (!valid_[i])
			return false;
	}
	return true;
}

Wt::WString BoundariesGroupBox::prepareLabel(std::string u, std::string v, std::string w)
{
	Wt::WString pattern("F({1}, {2}, {3}) =");
	return pattern.arg(u).arg(v).arg(w);
}

void BoundariesGroupBox::update(Limits limits)
{
	xLeftLabel_->setText(prepareLabel(toString(limits.xLeftLimit), "v", "w"));
	xRightLabel_->setText(prepareLabel(toString(limits.xRightLimit), "v", "w"));

	yLeftLabel_->setText(prepareLabel("u", toString(limits.yLeftLimit), "w"));
	yRightLabel_->setText(prepareLabel("u", toString(limits.yRightLimit), "w"));

	zLeftLabel_->setText(prepareLabel("u", "v", toString(limits.zLeftLimit)));
	zRightLabel_->setText(prepareLabel("u", "v", toString(limits.zRightLimit)));
}

void BoundariesGroupBox::enableAll()
{
	enabled_ = true;
	std::for_each(this->children().begin(), this->children().end(), [&](WWidget* w) { w->enable(); });
	if (areAllValid())
		validSignal_();
}

void BoundariesGroupBox::disableAll()
{
	enabled_ = false;
	std::for_each(this->children().begin(), this->children().end(), [&](WWidget* w) { w->disable(); });
}

std::array<PolynomialData, 6> BoundariesGroupBox::boundaries()
{
	std::array<PolynomialData, 6> result;
	result[0] = xLeft_->value();
	result[1] = xRight_->value();
	result[2] = yLeft_->value();
	result[3] = yRight_->value();
	result[4] = zLeft_->value();
	result[5] = zRight_->value();

	return result;
}

BoundariesGroupBox::~BoundariesGroupBox()
{
	my_delete(xLeft_);
	my_delete(xLeftLabel_);
	my_delete(xRight_);
	my_delete(xRightLabel_);

	my_delete(yLeft_);
	my_delete(yLeftLabel_);
	my_delete(yRight_);
	my_delete(yRightLabel_);

	my_delete(zLeft_);
	my_delete(zLeftLabel_);
	my_delete(zRight_);
	my_delete(zRightLabel_);
}

///// RESULT BOX

RhsGroupBox::RhsGroupBox(WContainerWidget* parent) : WGroupBox("Right hand side", parent)
{
	this->addStyleClass("child-legend");
	initComponents();
	initEvents();
}

void RhsGroupBox::initComponents()
{
	Wt::WTemplate* t = new Wt::WTemplate(Wt::WString::tr("rhs-form"), this);
	t->addFunction("id", &Wt::WTemplate::Functions::id);
	rhs_ = new PolynomialLineEdit();
	t->bindWidget("rhs", rhs_);
	rhs_->setDisabled(true);
}

void RhsGroupBox::initEvents()
{
	rhs_->blurred().connect(std::bind([=](){
		if (rhs_->validate() == Wt::WValidator::Valid)
			validSignal_();
		else
			invalidSignal_();
	}));
}

void RhsGroupBox::enable()
{
	rhs_->setDisabled(false);
	if (rhs_->validate() == Wt::WValidator::Valid)
		validSignal_();
}

RhsGroupBox::~RhsGroupBox()
{
	my_delete(rhs_);
}

///// RESULT BOX

SystemSettingsBox::SystemSettingsBox(WContainerWidget* parent) : WGroupBox("Chart settings", parent)
{
	setMargin(20, Wt::Bottom);
	addStyleClass("main-legend");

	initComponents();
}

void SystemSettingsBox::initComponents()
{
	density_ = new Wt::WComboBox();
	for (int i = 0; i < stepsArraySize; i++)
	{
		density_->addItem(std::to_string(possibleDensities[i]));
	}
	density_->setCurrentIndex(3);

	Wt::WTemplate* t = new Wt::WTemplate(Wt::WString::tr("density-form"), this);
	t->bindWidget("density", density_);
}

unsigned SystemSettingsBox::density()
{
	return possibleDensities[density_->currentIndex()];
}

double SystemSettingsBox::stepSize()
{
	return possibleStepSizes[density_->currentIndex()];
}

SystemSettingsBox::~SystemSettingsBox()
{
	my_delete(density_);
}