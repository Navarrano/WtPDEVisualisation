#include "InputPanel.h"

#include <Wt/WText>
#include <Wt/WGroupBox>
#include <Wt/WLabel>
#include <Wt/WTemplate>

std::string toString(const double val, short precision = 2)
{
	std::stringstream str;
	str << std::fixed << std::setprecision(precision) << val;
	return str.str();
}


InputPanel::InputPanel(WContainerWidget* parent) : AbstractPanel("Input data", parent)
{
	initComponents();
	initEvents();

	std::ifstream ifs("F:/LaplaceSoln2.dat");
	//ifs >> pdeData_;
}

void InputPanel::initComponents()
{
	WGroupBox* pdeBox = new WGroupBox("Partial derivative equation", root());
	pdeBox->addStyleClass("main-legend");

	WGroupBox* ordersBox = new WGroupBox("Orders of surfaces", pdeBox);
	ordersBox->setStyleClass("child-legend");
	WTemplate* t = new WTemplate(WString::tr("surfaces-form"), ordersBox);
	t->addFunction("id", &WTemplate::Functions::id);
	xOrder_ = initComboBox();
	t->bindWidget("xSurface", xOrder_);
	yOrder_ = initComboBox();
	t->bindWidget("ySurface", yOrder_);
	zOrder_ = initComboBox();
	t->bindWidget("zSurface", zOrder_);

	limitsBox_ = new LimitsGroupBox(pdeBox);
	boundariesBox_ = new BoundariesGroupBox(pdeBox);
	
	WGroupBox* resultBox = new WGroupBox("Right hand side", root());
	resultBox->setStyleClass("child-legend");
	t = new WTemplate(WString::tr("rhs-form"), resultBox);
	t->addFunction("id", &WTemplate::Functions::id);
	rhs_ = new PolynomialLineEdit();
	t->bindWidget("rhs", rhs_);

	WGroupBox *groupBox = new WGroupBox("System settings", root());
	groupBox->setMargin(20, Bottom);
	groupBox->addStyleClass("main-legend");
	density_ = new WComboBox();
	for (int i = 0; i < stepsArraySize; i++)
	{
		density_->addItem(std::to_string(possibleDensities[i]));
	}
	density_->setCurrentIndex(3);

	t = new WTemplate(WString::tr("density-form"), root());
	t->bindWidget("density", density_);

	new WText("<hr/>", root());
	calculateButton_ = new WPushButton("Calculate", root());
}

void InputPanel::initEvents()
{
	

	calculateButton_->clicked().connect(std::bind([=](){
		ChartData chartData;
		unsigned density = possibleDensities[density_->currentIndex()];
		pdeData_ = calculatePDE(SurfaceOrders(3, 3, 3), limitsBox_->getLimits(), boundariesBox_->boundaries(), rhs_->value());
		chartData.xStart = pdeData_.GetLeftLimitU();
		chartData.xEnd = pdeData_.GetRightLimitU();
		chartData.yStart = pdeData_.GetLeftLimitV();
		chartData.yEnd = pdeData_.GetRightLimitV();
		chartData.zStart = pdeData_.GetLeftLimitW();
		chartData.zEnd = pdeData_.GetRightLimitW();
		chartData.xSize = (chartData.xEnd - chartData.xStart) * density + 1;
		chartData.ySize = (chartData.yEnd - chartData.yStart) * density + 1;
		chartData.zSize = (chartData.zEnd - chartData.zStart) * density + 1;
		chartData.totalNbPts = chartData.xSize * chartData.ySize * chartData.zSize;
		chartData.points = pdeData_.ComputePoints(chartData.xSize - 1, chartData.ySize - 1, chartData.zSize - 1);
		chartData.step = possibleStepSizes[density_->currentIndex()];
		activate(chartData);
	}));

	limitsBox_->valid().connect(std::bind([=](){
		if (!boundariesBox_->enabled())
		{
			boundariesBox_->enableAll();
			rhs_->setDisabled(false);
			calculateButton_->setDisabled(false);
		}
		boundariesBox_->update(limitsBox_->getLimits());
	}));
	limitsBox_->invalid().connect(std::bind([=](){
		if (boundariesBox_->enabled())
		{
			boundariesBox_->disableAll();
			rhs_->setDisabled(true);
			calculateButton_->setDisabled(true);
		}
	}));
}

FBspVol& InputPanel::preperVolData(SurfaceOrders orders, Limits limits)
{
	short orderx = orders.xOrder;
	short ordery = orders.yOrder;
	short orderz = orders.zOrder;
	short segments = 5;

	double stepX = (limits.xRightLimit - limits.xLeftLimit) / (double)segments;
	double stepY = (limits.yRightLimit - limits.yLeftLimit) / (double)segments;
	double stepZ = (limits.zRightLimit - limits.zLeftLimit) / (double)segments;

	short numx = orderx + segments - 1;
	short numy = ordery + segments - 1;
	short numz = orderz + segments - 1;

	// knot sets
	Vector<double> knotX(numx + orderx), knotY(numy + ordery), knotZ(numz + orderz);

	for (int i = 0; i<orderx; i++) knotX[i] = limits.xLeftLimit;
	for (int i = orderx; i<numx; i++) knotX[i] = limits.xLeftLimit + (i - orderx + 1)*stepX;
	for (int i = numx; i<numx + orderx; i++) knotX[i] = limits.xRightLimit;

	for (int i = 0; i<ordery; i++) knotY[i] = limits.yLeftLimit;
	for (int i = ordery; i<numy; i++) knotY[i] = limits.yLeftLimit + (i - ordery + 1)*stepY;
	for (int i = numy; i<numy + ordery; i++) knotY[i] = limits.yRightLimit;

	for (int i = 0; i<orderz; i++) knotZ[i] = limits.zLeftLimit;
	for (int i = orderz; i<numz; i++) knotZ[i] = limits.zLeftLimit + (i - orderz + 1)*stepZ;
	for (int i = numz; i<numz + orderz; i++) knotZ[i] = limits.zRightLimit;

	Matrix3D<double> cpts(numx, numy, numz, 0.0);

	return FBspVol(cpts, knotX, knotY, knotZ, orderx, ordery, orderz, numx, numy, numz);
}

FBspVol& InputPanel::calculatePDE(SurfaceOrders surfaceOrders, Limits limits, std::array<PolynomialData, 6> boundaries, PolynomialData rhs)
{
	FBspVol b = preperVolData(surfaceOrders, limits);

	int num_pointforce = 0, num_distforce = 1, num_supports = 0;

	Vector<int> nums(40);

	nums[0] = num_pointforce;
	nums[1] = num_distforce;
	nums[2] = num_supports;

	Vector<double> pointforce(4 * num_pointforce);
	Vector<double> distforce(6 * num_distforce);
	Vector<double> supports(4 * num_supports);
	Vector<int> orders(3 * num_distforce);

	for (int i = 0; i < 4 * num_pointforce; i++) pointforce[i] = 0;
	distforce[0] = limits.xLeftLimit;
	distforce[1] = limits.xRightLimit;
	distforce[2] = limits.yLeftLimit;
	distforce[3] = limits.yRightLimit;
	distforce[4] = limits.zLeftLimit;
	distforce[5] = limits.zRightLimit;

	for (int i = 0; i < 3 * num_distforce; i++) orders[i] = rhs.order_;

	Vector<Matrix3D<double> > vol(num_distforce);
	int ind = 0;

	if (num_distforce > 0) {
		for (int i = 0; i<num_distforce; i++) {
			Matrix3D<double> voldata(orders[ind], orders[ind + 1], orders[ind + 2]);
			int l = 0;
			for (int k = 0; k<orders[ind + 2]; k++)
				for (int i = 0; i<orders[ind]; i++)
					for (int j = 0; j<orders[ind + 1]; j++) voldata[k][i][j] = rhs.coeffs_[l++];
			ind = ind + 3;
			vol[i] = voldata;
		}
	}

	for (int i = 0; i<4 * num_supports; i++) supports[i] = 0;
	double L = limits.xRightLimit - limits.xLeftLimit;
	double W = limits.yRightLimit - limits.yLeftLimit;
	double H = limits.zRightLimit - limits.zLeftLimit;

	// boundary conditions

	Vector<int> natgeom(6);
	for (int i = 0; i < 6; i++) natgeom[i] = 0;

	if (rhs.order_ == 1)
		return b.ComputeSolution(L, W, H, nums, supports, pointforce, distforce, orders, vol, natgeom, bound);
	else
		return b.ComputeSolution1(L, W, H, nums, supports, pointforce, distforce, orders, vol, natgeom, bound);
}

WComboBox* InputPanel::initComboBox()
{
	WComboBox* cb = new WComboBox();
	for (int i = 2; i <= 6; i++)
	{
		cb->addItem(std::to_string(i));
	}
	cb->setCurrentIndex(1);
	return cb;
}

InputPanel::~InputPanel()
{
	delete calculateButton_;
	delete limitsBox_;
	delete boundariesBox_;
	delete rhs_;
	delete density_;
	delete xOrder_;
	delete yOrder_;
	delete zOrder_;
}

// CUSTOM VALIDATORS 
CustomValidator::CustomValidator(DoubleLineEdit* compareWith, WObject* parent) : WDoubleValidator(parent), compareWith_(compareWith) {}

WValidator::Result CustomValidator::validate(const WString &input) const
{
	WValidator::Result r = WDoubleValidator::validate(input);
	if (r.state() != WValidator::Valid)
		return r;
	if (compareWith_->validate() != WValidator::Valid)
		return WValidator::Result(WValidator::Invalid, "Minimum value is required");
	double val = compareWith_->value();
	if (std::stod(input.toUTF8()) <= val)
		return WValidator::Result(WValidator::Invalid, "Value has to be greater than minimum value");
	else
		return WValidator::Result(WValidator::Valid);
}

PolynomialValidator::PolynomialValidator(WObject* parent) : WValidator(parent), twoVar_(false) {};
PolynomialValidator::PolynomialValidator(const char invalidVariable, WObject* parent) : WValidator(parent), twoVar_(true), invalidVariable_(invalidVariable) {};

WValidator::Result PolynomialValidator::validate(const WString &input) const
{
	WValidator::Result r = WValidator::validate(input);
	if (r.state() != WValidator::Valid)
		return r;
	PolynomialParser* p;
	if (twoVar_)
		p = new PolynomialParser(input.toUTF8(), invalidVariable_);
	else
		p =  new PolynomialParser(input.toUTF8());
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

// DOUBLE LINE EDIT

DoubleLineEdit::DoubleLineEdit(WContainerWidget* parent) : WLineEdit(parent)
{
	WValidator* validator = new WDoubleValidator(this);
	validator->setMandatory(true);
	this->setValidator(validator);
	WLineEdit::changed().connect(std::bind([=]() { signal_(this); }));
}

DoubleLineEdit::DoubleLineEdit(DoubleLineEdit* compareWith, WContainerWidget* parent)
{
	WValidator* validator = new CustomValidator(compareWith);
	validator->setMandatory(true);
	this->setValidator(validator);
	WLineEdit::blurred().connect(std::bind([=]() { signal_(this); }));
}

const double DoubleLineEdit::value()
{
	if (validate() == WValidator::Valid)
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

DoubleLineEdit::~DoubleLineEdit() {}

/////// POLYNOMIAL LINE EDIT

PolynomialLineEdit::PolynomialLineEdit(WContainerWidget* parent) : WLineEdit(parent), twoVar_(false)
{
	WValidator* validator = new PolynomialValidator(this);
	validator->setMandatory(true);
	this->setValidator(validator);
	WLineEdit::blurred().connect(std::bind([=]() { this->validate(); signal_(this); }));
}

PolynomialLineEdit::PolynomialLineEdit(const char invalidVariable, WContainerWidget* parent) : WLineEdit(parent), twoVar_(true), invalidVariable_(invalidVariable)
{
	WValidator* validator = new PolynomialValidator(invalidVariable_, this);
	validator->setMandatory(true);
	this->setValidator(validator);
	WLineEdit::blurred().connect(std::bind([=]() { this->validate(); signal_(this); }));
}
PolynomialData PolynomialLineEdit::value()
{
	PolynomialParser* p;
	if (twoVar_)
		p = new PolynomialParser((this->valueText()).toUTF8(), invalidVariable_);
	else
		p = new PolynomialParser((this->valueText()).toUTF8());
	return p->parse();
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
	WTemplate* t = new WTemplate(WString::tr("limits-form"), this);
	t->addFunction("id", &WTemplate::Functions::id);

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
			if (inputs_[i]->validate() == WValidator::Valid)
			{
				valid_[i] = true;
				if (i <= 2)
				{
					if (inputs_[i + 3]->validate() == WValidator::Valid)
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
}

////// BOUNDARIES GROUP BOX

BoundariesGroupBox::BoundariesGroupBox(WContainerWidget* parent) : WGroupBox("Boundaries", parent), enabled_(true)
{
	this->addStyleClass("child-legend");
	
	initComponents();
	initEvents();
}

void BoundariesGroupBox::initComponents()
{
	WTemplate* t = new WTemplate(WString::tr("boundaries-form"), this);
	t->addFunction("id", &WTemplate::Functions::id);

	xLeft_ = new PolynomialLineEdit('u');
	xLeftLabel_ = new WText("F(0.00 ,v, w) =");
	xRight_ = new PolynomialLineEdit('u');
	xRightLabel_ = new WText("F(1.00 ,v ,w) =");

	yLeft_ = new PolynomialLineEdit('v');
	yLeftLabel_ = new WText("F(u, 0.00, w) =");
	yRight_ = new PolynomialLineEdit('v');
	yRightLabel_ = new WText("F(u, 1.00, w) =");

	zLeft_ = new PolynomialLineEdit('w');
	zLeftLabel_ = new WText("F(u, v, 0.00) =");
	zRight_ = new PolynomialLineEdit('w');
	zRightLabel_ = new WText("F(u, v, 1.00) =");

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

}

WString BoundariesGroupBox::prepareLabel(std::string u, std::string v, std::string w)
{
	WString pattern("F({1}, {2}, {3}) =");
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
	delete xLeft_;
	delete xLeftLabel_;
	delete xRight_;
	delete xRightLabel_;

	delete yLeft_;
	delete yLeftLabel_;
	delete yRight_;
	delete yRightLabel_;

	delete zLeft_;
	delete zLeftLabel_;
	delete zRight_;
	delete zRightLabel_;
}



///// RESULT BOX

ResultBox::ResultBox(WContainerWidget* parent) : WGroupBox("Left hands side", parent)
{
	this->addStyleClass("child-legend");

	
}

ResultBox::~ResultBox()
{
	delete lhs_;
}