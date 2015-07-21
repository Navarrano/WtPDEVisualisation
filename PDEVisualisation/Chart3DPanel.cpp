#include "Chart3DPanel.h"

#include <Wt/WText>
#include <Wt/WHBoxLayout>
#include <Wt/WLabel>
#include <Wt/WTemplate>

#include "ViewButtonsWidget.h"


WTemplate* preparseSingleRowFormTemplate(WFormWidget* input, std::string labelText, int labelColSize = 2, int labelOffset = 2, int inputColSize = 4)
{
	WTemplate* t = new WTemplate(WString::tr("single-row-form"));
	t->bindWidget("input", input);
	t->bindString("label", labelText);
	t->bindInt("labelColSize", labelColSize);
	t->bindInt("labelOffset", labelOffset);
	t->bindInt("inputColSize", inputColSize);

	return t;
}


Chart3DPanel::Chart3DPanel(WContainerWidget *parent) : AbstractPanel("3D Volume Chart", parent), scatter_(0)
{
	this->collapse();

	initComponents();
	initEvents();
}

void Chart3DPanel::initComponents()
{
	WHBoxLayout *hbox = new WHBoxLayout();
	WContainerWidget *chartContainer = new WContainerWidget();
	WContainerWidget *settingsContainer = new WContainerWidget();

	volumeChart_ = new BaseChart(chartContainer);
	initChartLook();

	settingsBox_ = new SettingsBox(settingsContainer);

	ViewButtonsWidget *viewButtons = new ViewButtonsWidget(chartContainer);
	viewButtons->clicked().connect(boost::bind(&BaseChart::changeView, volumeChart_, _1));

	hbox->addWidget(chartContainer);
	hbox->addWidget(settingsContainer, 1);

	root()->setLayout(hbox);
}

void Chart3DPanel::initEvents()
{
	settingsBox_->changedRate().connect(std::bind([=](DisplayRate xRate, DisplayRate yRate, DisplayRate zRate){
		volumeData_->setXRate(xRate);
		volumeData_->setYRate(yRate);
		volumeData_->setZRate(zRate);
		volumeChart_->updateChart(GLContext | GLTextures);
	}, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));

	settingsBox_->changed().connect(std::bind([=](VolumeSettings settings){
		volumeData_->setXRate(settings.xRate);
		volumeData_->setYRate(settings.yRate);
		volumeData_->setZRate(settings.zRate);
		volumeData_->setPointSize(settings.generalSize);
		volumeData_->setSkippedPointSize(settings.skippedSize);
		volumeChart_->updateChart(GLContext | GLTextures);
	}, std::placeholders::_1));

	settingsBox_->clip().connect(std::bind([=](Limits limits){
		volumeData_->setLimits(limits);
		volumeChart_->updateChart(GLContext | GLTextures);
	}, std::placeholders::_1));

	settingsBox_->getCheckBox()->checked().connect(std::bind([=](){
		volumeChart_->toggleColorMap(true);
	}));
	settingsBox_->getCheckBox()->unChecked().connect(std::bind([=](){
		volumeChart_->toggleColorMap(false);
	}));
}

void Chart3DPanel::initChartLook()
{
	volumeChart_->resize(826, 826);
	volumeChart_->setBackground(WColor(210, 210, 210));
	WFont labelFont = volumeChart_->axis(Axis::XAxis_3D).labelFont();
	WFont titleFont = volumeChart_->axis(Axis::XAxis_3D).titleFont();
	titleFont.setWeight(WFont::Bold);
	labelFont.setSize(WFont::Medium);
	titleFont.setSize(WFont::Large);

	WPen textPen(WColor(232, 73, 12));

	volumeChart_->axis(Axis::XAxis_3D).setLabelFont(labelFont);
	volumeChart_->axis(Axis::XAxis_3D).setTitleFont(titleFont);
	volumeChart_->axis(Axis::XAxis_3D).setTextPen(textPen);

	volumeChart_->axis(Axis::YAxis_3D).setLabelFont(labelFont);
	volumeChart_->axis(Axis::YAxis_3D).setTitleFont(titleFont);
	volumeChart_->axis(Axis::YAxis_3D).setTextPen(textPen);

	volumeChart_->axis(Axis::ZAxis_3D).setLabelFont(labelFont);
	volumeChart_->axis(Axis::ZAxis_3D).setTitleFont(titleFont);
	volumeChart_->axis(Axis::ZAxis_3D).setTextPen(textPen);
}

void Chart3DPanel::update(ChartData& chartData)
{
	volumeData_ = new VolumeData(chartData, volumeChart_);

	volumeChart_->clearDatasets();
	scatter_ = volumeChart_->addVolumeDataset(volumeData_);

	settingsBox_->update(chartData);
	AbstractPanel::update(chartData);
}

Chart3DPanel::~Chart3DPanel()
{
	delete volumeChart_;
	delete settingsBox_;
}

////////////////// SETTINGS BOX ////////////////////////

SettingsBox::SettingsBox(WContainerWidget *parent) : WContainerWidget(parent)
{
	initComponenets();
	initEvents();
}

void SettingsBox::initComponenets()
{
	WGroupBox* radiosGroup = new WGroupBox("Points frequency", this);
	xRate_ = new DisplayRadioButtons("U", radiosGroup);
	yRate_ = new DisplayRadioButtons("V", radiosGroup);
	zRate_ = new DisplayRadioButtons("W", radiosGroup);

	clippingBox_ = new ClippingGroupBox(this);
	clippingBox_->setMargin(15, Top);

	WGroupBox* sizesGroup = new WGroupBox("Points sizes", this);
	sizesGroup->setMargin(15, Top);
	WLabel* generalLabel = new WLabel("General: ", sizesGroup);
	general_ = new WComboBox(sizesGroup);
	general_->setDisabled(true);
	initComboBox(general_, 5);
	generalLabel->setBuddy(general_);
	WLabel* skippedLabel = new WLabel("Skipped: ", sizesGroup);
	skipped_ = new WComboBox(sizesGroup);
	skipped_->setDisabled(true);
	initComboBox(skipped_, 0);
	skippedLabel->setBuddy(skipped_);

	new WText("<hr/>", this);

	showColorMap_ = new WCheckBox("Show color map", this);
	showColorMap_->setChecked(true);
}

void SettingsBox::initEvents()
{
	xRate_->changed().connect(std::bind([=]() {
		settingsChangedSignal_(generateSettings());
	}));
	yRate_->changed().connect(std::bind([=]() {
		settingsChangedSignal_(generateSettings());
	}));
	zRate_->changed().connect(std::bind([=]() {
		settingsChangedSignal_(generateSettings());
	}));
	general_->changed().connect(std::bind([=](){
		settingsChangedSignal_(generateSettings());
	}));
	skipped_->changed().connect(std::bind([=](){
		settingsChangedSignal_(generateSettings());
	}));
}

boost::signal<void(Limits)>& SettingsBox::clip()
{
	return clippingBox_->clicked();
}

VolumeSettings SettingsBox::generateSettings()
{
	VolumeSettings settings;
	settings.xRate = xRate_->getRate();
	settings.yRate = yRate_->getRate();
	settings.zRate = zRate_->getRate();
	settings.generalSize = general_->currentIndex() + 1;
	settings.skippedSize = skipped_->currentIndex() + 1;

	return settings;
}

void SettingsBox::initComboBox(WComboBox* cb, int index)
{
	for (int i = 1; i <= 15; i++)
	{
		cb->addItem(std::to_string(i));
	}
	cb->setCurrentIndex(index);
}

void SettingsBox::update(ChartData& chartData)
{
	activateControls();
	general_->setCurrentIndex(5);
	skipped_->setCurrentIndex(0);
	clippingBox_->update(chartData);
}

void SettingsBox::activateControls()
{
	xRate_->activate();
	yRate_->activate();
	zRate_->activate();
	skipped_->setDisabled(false);
	general_->setDisabled(false);
}

SettingsBox::~SettingsBox() 
{
	delete xRate_;
	delete yRate_;
	delete zRate_;
	delete general_;
	delete skipped_;
	delete clippingBox_;
	delete showColorMap_;
}

////////////////// DISPLAY RADIO BUTTONS ////////////////////////

DisplayRadioButtons::DisplayRadioButtons(char* label, WContainerWidget *parent) : WContainerWidget(parent)
{
	initComponents(label);
	initEvents();
}
DisplayRadioButtons::~DisplayRadioButtons()
{
	delete group_;
	delete everyPoint_;
	delete everySecondPoint_;
	delete everyFifthPoint_;
	delete everyTenthPoint_;
}

void DisplayRadioButtons::initComponents(char* text)
{
	char str[6];
	strcpy_s(str, text);
	strcat_s(str, ":  ");
	WText *label = new WText(str, this);

	group_ = new WButtonGroup(this);
	everyPoint_ = new WRadioButton("Every", this);
	group_->addButton(everyPoint_, 0);
	everySecondPoint_ = new WRadioButton("Every second", this);
	group_->addButton(everySecondPoint_, 1);
	everyFifthPoint_ = new WRadioButton("Every fifth", this);
	group_->addButton(everyFifthPoint_, 2);
	everyTenthPoint_ = new WRadioButton("Every tenth", this);
	group_->addButton(everyTenthPoint_, 3);
	group_->setSelectedButtonIndex(0);

	everyPoint_->setDisabled(true);
	everySecondPoint_->setDisabled(true);
	everyFifthPoint_->setDisabled(true);
}

void DisplayRadioButtons::initEvents()
{
	group_->checkedChanged().connect(std::bind([=](WRadioButton* selection){
		radioChangedSignal_();
	}, std::placeholders::_1));
}

DisplayRate DisplayRadioButtons::getRate()
{
	switch (group_->selectedButtonIndex())
	{
	case 0:
		return DisplayRate::Every;
		break;
	case 1:
		return DisplayRate::EverySecond;
		break;
	case 2:
		return DisplayRate::EveryFifth;
		break;
	case 3:
		return DisplayRate::EveryTenth;
		break;
	default:
		break;
	}
}

void DisplayRadioButtons::activate()
{
	everyPoint_->setDisabled(false);
	everySecondPoint_->setDisabled(false);
	everyFifthPoint_->setDisabled(false);
	everyTenthPoint_->setDisabled(false);

	group_->setSelectedButtonIndex(0);
}

//////// CLIPPING SLIDERS ////////

ClippingAxisSliders::ClippingAxisSliders(char* label, WContainerWidget* parent) : WContainerWidget(parent)
{
	initComponents();
	initEvents();
}

void ClippingAxisSliders::update(double min, double step, unsigned nbPts)
{
	min_ = min;
	step_ = step;

	leftLimitSlider_->setDisabled(false);
	leftLimitSlider_->setRange(0, nbPts);
	leftLimitSlider_->setTickInterval(5);
	leftLimitSlider_->setValue(0);
	leftLimitText_->setText(tickToTextValue(0));

	rightLimitSlider_->setDisabled(false);
	rightLimitSlider_->setRange(0, nbPts);
	rightLimitSlider_->setTickInterval(5);
	rightLimitSlider_->setValue(nbPts);
	rightLimitText_->setText(tickToTextValue(nbPts));
}

WString ClippingAxisSliders::tickToTextValue(int tick)
{
	std::ostringstream out;
	out << std::setprecision(2) << tick * step_ + min_;
	return WString(out.str());
}

double ClippingAxisSliders::getLeftLimit()
{
	return leftLimitSlider_->value() * step_ + min_;
}

double ClippingAxisSliders::getRightLimit()
{
	return rightLimitSlider_->value() * step_ + min_;
}

void ClippingAxisSliders::initComponents()
{
	leftLimitSlider_ = initSlider();
	leftLimitText_ = new WText(this);
	leftLimitText_->setStyleClass("slider-text");
	leftLimitSlider_->setInline(false);
	rightLimitSlider_ = initSlider();
	rightLimitText_ = new WText(this);
	rightLimitText_->setStyleClass("slider-text");
	rightLimitSlider_->setInline(false);
}

void ClippingAxisSliders::initEvents()
{
	leftLimitSlider_->valueChanged().connect(std::bind([=](){
		if (leftLimitSlider_->value() >= rightLimitSlider_->value())
		{
			leftLimitSlider_->setValue(rightLimitSlider_->value());
		}
		leftLimitText_->setText(tickToTextValue(leftLimitSlider_->value()));
	}));

	rightLimitSlider_->valueChanged().connect(std::bind([=](){
		if (leftLimitSlider_->value() >= rightLimitSlider_->value())
		{
			rightLimitSlider_->setValue(leftLimitSlider_->value());
		}
		rightLimitText_->setText(tickToTextValue(rightLimitSlider_->value()));
	}));
}

WSlider* ClippingAxisSliders::initSlider()
{
	WSlider* slider = new WSlider(this);
	slider->resize(345, 50);
	slider->setTickPosition(WSlider::TicksAbove);
	slider->setDisabled(true);

	return slider;
}

ClippingAxisSliders::~ClippingAxisSliders()
{
	delete leftLimitSlider_;
	delete leftLimitText_;
	delete rightLimitSlider_;
	delete rightLimitText_;
}

//////// CLIPPING GROUP BOX ///////////

ClippingGroupBox::ClippingGroupBox(WContainerWidget* parent) : WGroupBox("Clipping", parent)
{
	initComponents();
	initEvents();
}

void ClippingGroupBox::initComponents()
{
	uSliders_ = new ClippingAxisSliders("U", this);
	new WText("<hr/>", this);
	vSliders_ = new ClippingAxisSliders("V", this);
	new WText("<hr/>", this);
	wSliders_ = new ClippingAxisSliders("W", this);
	new WText("<hr/>", this);

	clip_ = new WPushButton("Clip", this);
	clip_->setDisabled(true);
}

void ClippingGroupBox::initEvents()
{
	clip_->clicked().connect(std::bind([=](){
		Limits l;
		l.xLeftLimit = uSliders_->getLeftLimit();
		l.xRightLimit = uSliders_->getRightLimit();
		l.yLeftLimit = vSliders_->getLeftLimit();
		l.yRightLimit = vSliders_->getRightLimit();
		l.zLeftLimit = wSliders_->getLeftLimit();
		l.zRightLimit = wSliders_->getRightLimit();
		
		clipButtonSignal_(l);
	}));
}

void ClippingGroupBox::update(ChartData& chartData)
{
	uSliders_->update(chartData.xStart, chartData.step, chartData.xSize - 1);
	vSliders_->update(chartData.yStart, chartData.step, chartData.ySize - 1);
	wSliders_->update(chartData.zStart, chartData.step, chartData.zSize - 1);

	clip_->setDisabled(false);
}

ClippingGroupBox::~ClippingGroupBox()
{
	delete uSliders_;
	delete vSliders_;
	delete wSliders_;
	delete clip_;
}