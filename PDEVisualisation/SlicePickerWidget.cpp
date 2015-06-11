#include "SlicePickerWidget.h"

#include <Wt/WPushButton>
#include <Wt/WHBoxLayout>

SlicePickerWidget::SlicePickerWidget(WContainerWidget* parent) : WContainerWidget(parent)
{
	this->setStyleClass("form-group");

	comboBox_ = new WComboBox(this);
	comboBox_->addItem("u");
	comboBox_->addItem("v");
	comboBox_->addItem("w");
	comboBox_->setMargin(10, Right);
	comboBox_->setDisabled(true);

	spinBox_ = new WDoubleSpinBox(this);
	spinBox_->setMargin(10, Right);
	spinBox_->setDisabled(true);

	showButton_ = new WPushButton("Show", this);
	showButton_->setDisabled(true);
	errorText_ = new WText("", this);
	errorText_->addStyleClass("help-block");

	showButton_->clicked().connect(std::bind([=]() {
		if (spinBox_->validate() == WValidator::Valid)
		{
			errorText_->setText("");
			buttonClickedSignal((Dimension)comboBox_->currentIndex(), spinBox_->value());
		}
		else
		{
			errorText_->setText(Wt::WString::fromUTF8("Invalid spin box value!"));
		}	
	}));

	comboBox_->changed().connect(this, &SlicePickerWidget::updateSpinBox);
}

void SlicePickerWidget::updateSpinBox()
{
	switch (comboBox_->currentIndex())
	{
	case 0:
		spinBox_->setRange(uMin_, uMax_);
		if (spinBox_->validate() == WValidator::Invalid)
		{
			spinBox_->setValue(uMin_);
		}
		break;
	case 1:
		spinBox_->setRange(vMin_, vMax_);
		if (spinBox_->validate() == WValidator::Invalid)
		{
			spinBox_->setValue(vMin_);
		}
		break;
	case 2:
		spinBox_->setRange(wMin_, wMax_);
		if (spinBox_->validate() == WValidator::Invalid)
		{
			spinBox_->setValue(wMin_);
		}
		break;
	default:
		break;
	}
}

boost::signal<void(Dimension, double)>& SlicePickerWidget::clicked()
{
	return buttonClickedSignal;
}

void SlicePickerWidget::update(ChartData chartData, int precision)
{
	uMin_ = chartData.xStart;
	uMax_ = chartData.xEnd;
	vMin_ = chartData.yStart;
	vMax_ = chartData.yEnd;
	wMin_ = chartData.zStart;
	wMax_ = chartData.zEnd;

	comboBox_->setDisabled(false);
	spinBox_->setDisabled(false);
	showButton_->setDisabled(false);

	spinBox_->setSingleStep(chartData.step);
	spinBox_->setDecimals(precision);
	spinBox_->setRange(uMin_, uMax_);
	spinBox_->setValue(uMin_);
	comboBox_->setCurrentIndex(0);
}

SlicePickerWidget::~SlicePickerWidget()
{
	delete comboBox_;
	delete spinBox_;
	delete showButton_;
	delete errorText_;
}
